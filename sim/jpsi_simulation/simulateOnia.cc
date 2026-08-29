#include <array>
#include <string>
#include <numeric>
#include <map>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TSpline.h>
#include <TNtuple.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/HIInfo.h"

namespace
{
    enum decayer
    {
        kPythia8 = 0,
        kEvtGen
    };

    enum tunes
    {
        kMonash = 0,
        kCRMode0,
        kCRMode2,
        kCRMode3,
        kHeavyIon
    };

    enum processes
    {
        kSoftQCD = 0,
        kHardQCD
    };

    enum systems
    {
        kPP = 0,
        kPbPb,
        kOO,
        kNeNe
    };

    std::map<int, int> arrSystems{{kPP, 2212}, {kPbPb, 1000822080}, {kOO, 1000080160}, {kNeNe, 1000100200}};
    const std::array<int, 4> absPdgOpenCharm = {411, 421, 431, 4122};
}


template<typename T, typename PPythia>
bool isFromBeauty(T &mothers, PPythia& pythia)
{
    for(auto const& mom : mothers) {
      int absPdgMom = std::abs(pythia.event[mom].id());
      if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
        (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
        (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5) {
          return true;
      }
    }

  return false;
}

struct EventClassifiers
{
    double spherocity{-1.};
    double rt{-1.};
    double flattenicity{-1.};
    int nchMid{0};
    int nchTrans{0};
    int nchForward{0};
};

//__________________________________________________________________________________________________
template<typename PPythia>
EventClassifiers computeEventClassifiers(PPythia& pythia)
{
    EventClassifiers classifiers;

    double sumPtMid{0.};
    double minSpherocityNumerator{1.e30};
    double leadingPt{-1.};
    double leadingPhi{0.};

    constexpr int nPhiSteps = 180;
    constexpr int nFlattenicityPhiBins = 8;
    constexpr int nFlattenicityEtaRegions = 2;
    constexpr int nFlattenicityCells = nFlattenicityPhiBins * nFlattenicityEtaRegions;
    std::array<double, nFlattenicityCells> flattenicityCells{};

    for (auto iPart{2}; iPart < pythia.event.size(); ++iPart) {
        const auto& particle = pythia.event[iPart];
        if (!particle.isFinal() || !particle.isCharged()) {
            continue;
        }

        const double eta = particle.eta();
        const double pt = particle.pT();
        const double phi = particle.phi();

        // Mid-rapidity charged-particle selection for spherocity and leading-particle axis.
        if (std::abs(eta) < 0.8 && pt > 0.15) {
            ++classifiers.nchMid;
            sumPtMid += pt;

            if (pt > leadingPt) {
                leadingPt = pt;
                leadingPhi = phi;
            }
        }

        // Forward charged-particle map for flattenicity, using V0-like eta regions.
        int etaRegion = -1;
        if (eta > 2.8 && eta < 5.1) {
            etaRegion = 0;
        } else if (eta < -1.7 && eta > -3.7) {
            etaRegion = 1;
        }

        if (etaRegion >= 0) {
            ++classifiers.nchForward;
            double phiPositive = phi;
            while (phiPositive < 0.) {
                phiPositive += 2. * TMath::Pi();
            }
            while (phiPositive >= 2. * TMath::Pi()) {
                phiPositive -= 2. * TMath::Pi();
            }
            int phiBin = static_cast<int>(phiPositive / (2. * TMath::Pi() / nFlattenicityPhiBins));
            if (phiBin >= nFlattenicityPhiBins) {
                phiBin = nFlattenicityPhiBins - 1;
            }
            flattenicityCells[etaRegion * nFlattenicityPhiBins + phiBin] += 1.;
        }
    }

    // Transverse spherocity: S0 = pi^2/4 * min_n (sum_i |pT_i x n| / sum_i pT_i)^2.
    if (classifiers.nchMid >= 3 && sumPtMid > 0.) {
        for (int iStep = 0; iStep < nPhiSteps; ++iStep) {
            const double phiAxis = TMath::Pi() * static_cast<double>(iStep) / static_cast<double>(nPhiSteps);
            const double sinAxis = std::sin(phiAxis);
            const double cosAxis = std::cos(phiAxis);
            double numerator{0.};

            for (auto iPart{2}; iPart < pythia.event.size(); ++iPart) {
                const auto& particle = pythia.event[iPart];
                if (!particle.isFinal() || !particle.isCharged()) {
                    continue;
                }
                if (std::abs(particle.eta()) >= 0.8 || particle.pT() <= 0.15) {
                    continue;
                }
                numerator += std::abs(particle.px() * sinAxis - particle.py() * cosAxis);
            }
            if (numerator < minSpherocityNumerator) {
                minSpherocityNumerator = numerator;
            }
        }
        classifiers.spherocity = (TMath::Pi() * TMath::Pi() / 4.) * std::pow(minSpherocityNumerator / sumPtMid, 2.);
    }

    // RT numerator: charged multiplicity in the transverse region relative to the leading charged particle.
    // In post-processing, use RT = nchTrans / <nchTrans> for the chosen event sample.
    if (leadingPt > 0.) {
        for (auto iPart{2}; iPart < pythia.event.size(); ++iPart) {
            const auto& particle = pythia.event[iPart];
            if (!particle.isFinal() || !particle.isCharged()) {
                continue;
            }
            if (std::abs(particle.eta()) >= 0.8 || particle.pT() <= 0.15) {
                continue;
            }

            double deltaPhi = std::abs(particle.phi() - leadingPhi);
            while (deltaPhi > TMath::Pi()) {
                deltaPhi = std::abs(deltaPhi - 2. * TMath::Pi());
            }
            if (deltaPhi > TMath::Pi() / 3. && deltaPhi < 2. * TMath::Pi() / 3.) {
                ++classifiers.nchTrans;
            }
        }
        classifiers.rt = static_cast<double>(classifiers.nchTrans);
    }

    // Flattenicity = sigma(cell activity) / mean(cell activity) in V0-like cells.
    double meanActivity{0.};
    for (const auto& activity : flattenicityCells) {
        meanActivity += activity;
    }
    meanActivity /= static_cast<double>(nFlattenicityCells);

    if (meanActivity > 0.) {
        double variance{0.};
        for (const auto& activity : flattenicityCells) {
            variance += std::pow(activity - meanActivity, 2.);
        }
        variance /= static_cast<double>(nFlattenicityCells * nFlattenicityCells);
        classifiers.flattenicity = std::sqrt(variance) / meanActivity;
    }

    return classifiers;
}

//__________________________________________________________________________________________________
template<typename PPythia>
//void runSimulation(PPythia& pythia, TNtuple* tuplePairs, TH1D* histEvents, TH2D* histImpParvsNpart, std::map<int, TH3F*> histPtVsY, int nEvents)
void runSimulation(PPythia& pythia, int tune, TNtuple* tuplePairs, TH1D* histEvents, TH2D* histImpParvsNpart, int nEvents)
{
    //__________________________________________________________
    // perform the simulation
    for (auto iEvent{0}; iEvent < nEvents; ++iEvent)
    {
        if (iEvent > 0 && iEvent % 10000 == 0)
        {

            std::cout << "Events done : " << iEvent << std::endl;
        }
        if (!pythia.next())
        {
            continue;
        }

        histEvents->Fill(1.5f);
        const auto classifiers = computeEventClassifiers(pythia);
        double impactParameter, nPart;
        if (tune == kHeavyIon) {
            impactParameter = pythia.info.hiInfo->b();
            nPart = pythia.info.hiInfo->nPartTarg() + pythia.info.hiInfo->nPartProj();
            histImpParvsNpart -> Fill(impactParameter, nPart);
        }

        std::vector<Pythia8::Particle> jPsis;
        std::vector<bool> jPsiFromB{};
        std::vector<std::vector<int>> jPsiMothers{};

        double ptOnia, yOnia, etaOnia, phiOnia;

        int nFound{0};
        for(auto iPart{2}; iPart<pythia.event.size(); ++iPart)
        {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            if (absPdg == 443 || absPdg == 100443 || absPdg == 553 || absPdg == 100553 || absPdg == 200553) {
                auto mothers = pythia.event[iPart].motherList();
                size_t initMomSize{0u};
                while (mothers.size() > initMomSize) {
                    auto mothersOld = mothers;
                    initMomSize = mothersOld.size();
                    for (auto const& mom : mothersOld) {
                        auto mothersNew = pythia.event[mom].motherList();
                        for (auto const& momNew : mothersNew) {
                            if (std::find(mothers.begin(), mothers.end(), momNew) == mothers.end() && std::abs(pythia.event[momNew].id()) != 2212) {
                                mothers.push_back(momNew);
                            }
                        }
                    }
                }
                bool fromB = isFromBeauty(mothers, pythia);

                ptOnia = pythia.event[iPart].pT();
                yOnia = pythia.event[iPart].y();
                etaOnia = pythia.event[iPart].eta();
                phiOnia = pythia.event[iPart].phi();

                if (tune == kHeavyIon) {
                    tuplePairs->Fill(ptOnia, yOnia, etaOnia, phiOnia, fromB, absPdg, nPart,
                                     classifiers.spherocity, classifiers.rt, classifiers.flattenicity,
                                     classifiers.nchMid, classifiers.nchTrans, classifiers.nchForward);
                } else {
                    tuplePairs->Fill(ptOnia, yOnia, etaOnia, phiOnia, fromB, absPdg, 2,
                                     classifiers.spherocity, classifiers.rt, classifiers.flattenicity,
                                     classifiers.nchMid, classifiers.nchTrans, classifiers.nchForward);
                }
            }
        }
    }
}

//__________________________________________________________________________________________________
void simulateOnia(int nEvents, int tune, int process, bool usePtHardBins, float energy, int system, int seed, std::string outFileNameRoot)
{
    //__________________________________________________________
    // create and configure pythia generator

    Pythia8::Pythia pythia;
    if(process == kSoftQCD)
    {
        pythia.readString("SoftQCD:inelastic = on"); // we only simulate inelastic processes
    }
    else if(process == kHardQCD)
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
        pythia.readString("HardQCD:gg2ccbar = on");
        pythia.readString("HardQCD:gg2bbbar = on");
        pythia.readString("HardQCD:qqbar2ccbar = on");
        pythia.readString("HardQCD:qqbar2bbbar = on");
    }

    // tune for charmonia
    pythia.readString("CharmoniumShower:all = on");

    // set tune
    if(tune == kMonash)
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kHeavyIon) 
    {
        pythia.readString(Form("HeavyIon:mode = 1"));
        // old config
        //pythia.readString("HeavyIon:SigFitNGen = 0");
        //pythia.readString("HeavyIon:SigFitDefAvNDb = 0.92");
        //pythia.readString("HeavyIon:SigFitDefPar = 1.34,2.61,1.18");
        // config stored in https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/common/pythia8/generator/pythia8_OO_536.cfg
        pythia.readString("Beams:frameType = 1");
        pythia.readString("ParticleDecays:limitTau0 = on");
        pythia.readString("ParticleDecays:tau0Max = 10.");
        pythia.readString("HeavyIon:SigFitNGen = 0");
        pythia.readString("HeavyIon:SigFitDefPar = 2.15,18.42,0.33");
    } else {
        std::cout << "[ERROR] No tune configured" << std::endl;
        return;
    }

    if (usePtHardBins)
    {
        pythia.readString("PhaseSpace:pTHatMin = 20");
        pythia.readString("PhaseSpace:pTHatMax = 200");
    }

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.settings.mode("Beams:idA", arrSystems[system]);
    pythia.settings.mode("Beams:idB", arrSystems[system]);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    //__________________________________________________________
    // define outputs
    auto tuplePairs = new TNtuple("tuplePairs", "tuplePairs", "pTOnia:yOnia:etaOnia:phiOnia:fromB:absPdg:nPart:spherocity:rt:flattenicity:nchMid:nchTrans:nchForward");
    auto histNumEvents = new TH1D("histNumEvents", "", 2, 0., 2.);
    histNumEvents->GetXaxis()->SetBinLabel(1, "num events set");
    histNumEvents->GetXaxis()->SetBinLabel(2, "num events succeded");
    auto histImpParvsNpart = new TH2D("histImpParvsNpart", ";b (fm);<N_{Part}>", 200, 0, 20, 500, 0, 50);
    auto histXsec = new TH1D(Form("histXsec_%d", seed), ";cross section (mb)", 1, 0., 1.);
    histXsec->GetXaxis()->SetBinLabel(1, "");
    histNumEvents->SetBinContent(1, nEvents);
    //std::map<int, TH3F*> histPtVsY;
    //histPtVsY[443] = new TH3F("histPtVsY_443", ";#it{p}_{T} (GeV/#it{c}); #it{y}; fromB", 500, 0., 50., 200, -10., 10., 2, -0.5, 1.5);
    //histPtVsY[100443] = new TH3F("histPtVsY_100443", ";#it{p}_{T} (GeV/#it{c}); #it{y}; fromB", 500, 0., 50., 200, -10., 10., 2, -0.5, 1.5);
    runSimulation(pythia, tune, tuplePairs, histNumEvents, histImpParvsNpart, nEvents);
    histXsec->SetBinContent(1, pythia.info.sigmaGen());
    histXsec->SetBinError(1, pythia.info.sigmaErr());

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    tuplePairs->Write();
    histNumEvents->Write();
    histImpParvsNpart->Write();
    histXsec->Write();
    //histPtVsY[443]->Write();
    //histPtVsY[100443]->Write();
    outFile.Close();
}

//__________________________________________________________________________________________________
int main(int argc, char* argv[])
{
    std::cout << "Number of inputs " << argc << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << i << "/" << argc - 1 << " `" << argv[i] << "`" << std::endl;
    }

    // Input order:
    //   1 -> number of events
    //   2 -> tune: 0 Monash, 1 CRMode0, 2 CRMode2, 3 CRMode3, 4 HeavyIon
    //   3 -> process: 0 SoftQCD, 1 HardQCD
    //   4 -> use pT-hard bins: 0 false, 1 true
    //   5 -> collision energy in GeV
    //   6 -> system: 0 pp, 1 PbPb, 2 OO, 3 NeNe
    //   7 -> random seed
    //   8 -> output ROOT file name
    //   9 -> test flag, prints configuration and exits before running

    int nEvents = 100000000;
    int tune = kMonash;
    int process = kSoftQCD;
    bool usePtHardBins = false;
    float energy = 13000.f;
    int system = kPP;
    int seed = 0;
    std::string outFileNameRoot = "onia_output.root";

    if (argc >= 2) {
        nEvents = std::atoi(argv[1]);
    }
    if (argc >= 3) {
        tune = std::atoi(argv[2]);
    }
    if (argc >= 4) {
        process = std::atoi(argv[3]);
    }
    if (argc >= 5) {
        usePtHardBins = (std::atoi(argv[4]) != 0);
    }
    if (argc >= 6) {
        energy = std::atof(argv[5]);
    }
    if (argc >= 7) {
        system = std::atoi(argv[6]);
    }
    if (argc >= 8) {
        seed = std::atoi(argv[7]);
    }
    if (argc >= 9) {
        outFileNameRoot = argv[8];
    }

    std::cout << "Running simulateOnia with:" << std::endl;
    std::cout << "  nEvents        = " << nEvents << std::endl;
    std::cout << "  tune           = " << tune << std::endl;
    std::cout << "  process        = " << process << std::endl;
    std::cout << "  usePtHardBins  = " << usePtHardBins << std::endl;
    std::cout << "  energy         = " << energy << " GeV" << std::endl;
    std::cout << "  system         = " << system << std::endl;
    std::cout << "  seed           = " << seed << std::endl;
    std::cout << "  output file    = " << outFileNameRoot << std::endl;

    if (argc >= 10) {
        std::cout << "Running test only. Exiting before event generation." << std::endl;
        return 0;
    }

    simulateOnia(nEvents, tune, process, usePtHardBins, energy, system, seed, outFileNameRoot);

    std::cout << "*************************************************************************************************************" << std::endl;
    std::cout << "*                             EVENTS ARE GENERATED SUCCESSFULLY                                             *" << std::endl;
    std::cout << "*************************************************************************************************************" << std::endl;

    return 0;
}
