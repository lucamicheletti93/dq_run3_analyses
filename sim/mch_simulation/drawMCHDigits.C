// ============================================================
// drawMCHDigits.C
//
// Legge i digit MCH da mchdigits.root e crea una mappa ADC
// per ogni Detection Element (DE), usando la segmentation
// MCH di O2.
//
// Inoltre carica la geometria TGeo da o2sim_geometry.root.
//
// OUTPUT:
//   mch_pad_maps.root
//
// Nel file di output vengono salvati SOLO i TH2Poly.
//
// Nessun canvas viene creato e nessun Draw() viene eseguito.
//
// ESECUZIONE:
//
//   root
//   .x drawMCHDigits.C
//
// oppure:
//
//   .x drawMCHDigits.C("mchdigits.root","mch_pad_maps.root")
//
// IMPORTANTE:
//   NON usare .L drawMCHDigits.C++
//   perché nel tuo ambiente ACLiC non linka automaticamente
//   le librerie O2 MCH Mapping.
// ============================================================

#include <iostream>
#include <map>
#include <string>
#include <algorithm>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH2Poly.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"

// O2 MCH
#include "MCHMappingInterface/Segmentation.h"
#include "DetectorsBase/GeometryManager.h"


// ============================================================
// O2 installation
// ============================================================

const char* O2MCHLIB =
    "/Users/lucamicheletti/alice/sw/osx_arm64/"
    "O2/daily-20260224-0000-local1/lib/";


// ============================================================
// Main function
// ============================================================

void drawMCHDigits(
    const char* inputFile = "sim_rej_list_high_stat/mchdigits_filtered.root",
    const char* outputFile = "sim_rej_list_high_stat/mch_pad_maps_filtered.root",
    const char* geometryFile = "sim_rej_list_high_stat/o2sim_geometry.root")
{
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "                    MCH DIGIT MAP\n";
    std::cout << "============================================================\n";
    std::cout << "\n";


    // ========================================================
    // 1. Load MCH Mapping libraries
    // ========================================================

    std::string lib3 =
        std::string(O2MCHLIB) +
        "libO2MCHMappingImpl3.dylib";

    std::string lib4 =
        std::string(O2MCHLIB) +
        "libO2MCHMappingImpl4.dylib";

    std::cout << "Loading MCH Mapping libraries...\n";

    int ret3 = gSystem->Load(lib3.c_str());
    int ret4 = gSystem->Load(lib4.c_str());

    std::cout << "  Impl3: " << ret3 << "\n";
    std::cout << "  Impl4: " << ret4 << "\n";

    // In ROOT:
    //   0 = library loaded successfully
    //   1 = library was already loaded
    //  <0 = error
    //
    // Therefore 0 or 1 are both OK.

    if (ret3 < 0 || ret4 < 0) {

        std::cerr
            << "ERROR: cannot load MCH Mapping libraries.\n";

        return;
    }


    // ========================================================
    // 2. Load O2 geometry
    // ========================================================

    std::cout << "\n";
    std::cout << "Loading O2 geometry...\n";
    std::cout << "  File: "
              << geometryFile
              << "\n";

    if (!gGeoManager) {

        o2::base::GeometryManager::loadGeometry(
            geometryFile);
    }

    if (!gGeoManager) {

        std::cerr
            << "ERROR: TGeo geometry was not loaded.\n"
            << "       Geometry file: "
            << geometryFile
            << "\n";

        return;
    }

    std::cout
        << "Geometry loaded successfully.\n";

    std::cout
        << "Geometry manager: "
        << gGeoManager->GetName()
        << "\n";

    if (gGeoManager->GetTopVolume()) {

        std::cout
            << "Top volume: "
            << gGeoManager->GetTopVolume()->GetName()
            << "\n";
    }

    if (gGeoManager->GetListOfVolumes()) {

        std::cout
            << "Number of TGeo volumes: "
            << gGeoManager->GetListOfVolumes()->GetEntries()
            << "\n";
    }


    // ========================================================
    // 3. Open input ROOT file
    // ========================================================

    TFile* fin =
        TFile::Open(inputFile, "READ");

    if (!fin || fin->IsZombie()) {

        std::cerr
            << "ERROR: cannot open input file: "
            << inputFile
            << "\n";

        return;
    }

    std::cout
        << "\nInput file: "
        << inputFile
        << "\n";


    // ========================================================
    // 4. Get TTree
    // ========================================================

    TTree* tree = nullptr;

    fin->GetObject("o2sim", tree);

    if (!tree) {

        std::cerr
            << "ERROR: TTree 'o2sim' not found.\n";

        fin->Close();
        delete fin;

        return;
    }

    std::cout
        << "Tree found: "
        << tree->GetName()
        << "\n";

    std::cout
        << "Entries: "
        << tree->GetEntries()
        << "\n";


    if (tree->GetEntries() == 0) {

        std::cerr
            << "ERROR: tree contains no entries.\n";

        fin->Close();
        delete fin;

        return;
    }


    // ========================================================
    // 5. TTreeReader
    //
    // Branches:
    //
    // MCHDigit.mDetID : Int_t[MCHDigit_]
    // MCHDigit.mPadID : Int_t[MCHDigit_]
    // MCHDigit.mADC   : UInt_t[MCHDigit_]
    // ========================================================

    TTreeReader reader(tree);

    TTreeReaderArray<Int_t> detID(
        reader,
        "MCHDigit.mDetID");

    TTreeReaderArray<Int_t> padID(
        reader,
        "MCHDigit.mPadID");

    TTreeReaderArray<UInt_t> adc(
        reader,
        "MCHDigit.mADC");


    // ========================================================
    // 6. Read first entry
    // ========================================================

    if (!reader.Next()) {

        std::cerr
            << "ERROR: cannot read first tree entry.\n";

        fin->Close();
        delete fin;

        return;
    }


    const Long64_t nDigits =
        detID.GetSize();


    std::cout
        << "\nNumber of digits = "
        << nDigits
        << "\n";


    if (padID.GetSize() != nDigits ||
        adc.GetSize() != nDigits) {

        std::cerr
            << "ERROR: inconsistent array sizes:\n";

        std::cerr
            << "  detID = "
            << detID.GetSize()
            << "\n";

        std::cerr
            << "  padID = "
            << padID.GetSize()
            << "\n";

        std::cerr
            << "  adc   = "
            << adc.GetSize()
            << "\n";

        fin->Close();
        delete fin;

        return;
    }


    // ========================================================
    // 7. Accumulate ADC per DE and pad
    //
    // padADC[DE][padID] = total ADC
    // ========================================================

    std::map<int, std::map<int, double>> padADC;


    std::cout
        << "\nReading digits...\n";


    for (Long64_t i = 0; i < nDigits; ++i) {

        const int deid =
            detID[i];

        const int pid =
            padID[i];

        const double charge =
            static_cast<double>(adc[i]);

        padADC[deid][pid] += charge;
    }


    // ========================================================
    // 8. Print DE summary
    // ========================================================

    std::cout
        << "\nDetection Elements found: "
        << padADC.size()
        << "\n";


    for (const auto& deEntry : padADC) {

        const int deid =
            deEntry.first;

        const auto& pads =
            deEntry.second;

        std::cout
            << "  DE "
            << deid
            << " : "
            << pads.size()
            << " active pads\n";
    }


    // ========================================================
    // 9. Open output ROOT file
    // ========================================================

    TFile* fout =
        TFile::Open(outputFile, "RECREATE");


    if (!fout || fout->IsZombie()) {

        std::cerr
            << "ERROR: cannot create output file: "
            << outputFile
            << "\n";

        fin->Close();
        delete fin;

        return;
    }


    // ========================================================
    // 10. Create one TH2Poly per DE
    //
    // NOTE:
    // At this stage the bins are still rectangles generated
    // from the O2 segmentation.
    //
    // The TGeo geometry has been loaded and is available
    // through gGeoManager.
    //
    // In the next step we can replace these rectangles with
    // the REAL pad polygons extracted from TGeo.
    // ========================================================

    int nDE = 0;
    int nInvalidPads = 0;


    for (const auto& deEntry : padADC) {

        const int deid =
            deEntry.first;

        const auto& pads =
            deEntry.second;


        std::cout << "\n";
        std::cout
            << "------------------------------------------------------------\n";

        std::cout
            << "DE "
            << deid
            << "\n";

        std::cout
            << "------------------------------------------------------------\n";


        // ====================================================
        // Get O2 segmentation
        // ====================================================

        const o2::mch::mapping::Segmentation& seg =
            o2::mch::mapping::segmentation(deid);


        std::cout
            << "Total pads in segmentation: "
            << seg.nofPads()
            << "\n";


        // ====================================================
        // Determine geometrical limits
        // ====================================================

        double xmin =  1.e30;
        double xmax = -1.e30;

        double ymin =  1.e30;
        double ymax = -1.e30;


        int nValidPads = 0;


        for (const auto& padEntry : pads) {

            const int pid =
                padEntry.first;


            if (!seg.isValid(pid)) {

                ++nInvalidPads;

                continue;
            }


            const double x =
                seg.padPositionX(pid);

            const double y =
                seg.padPositionY(pid);

            const double dx =
                seg.padSizeX(pid) / 2.0;

            const double dy =
                seg.padSizeY(pid) / 2.0;


            xmin =
                std::min(xmin, x - dx);

            xmax =
                std::max(xmax, x + dx);

            ymin =
                std::min(ymin, y - dy);

            ymax =
                std::max(ymax, y + dy);


            ++nValidPads;
        }


        if (nValidPads == 0) {

            std::cerr
                << "WARNING: DE "
                << deid
                << " contains no valid pads.\n";

            continue;
        }


        // ====================================================
        // Add small margin
        // ====================================================

        const double marginX =
            0.02 * (xmax - xmin);

        const double marginY =
            0.02 * (ymax - ymin);


        xmin -= marginX;
        xmax += marginX;

        ymin -= marginY;
        ymax += marginY;


        // ====================================================
        // Create TH2Poly
        // ====================================================

        TString hname;

        hname.Form(
            "hMCHADC_DE%d",
            deid);


        TString htitle;

        htitle.Form(
            "MCH DE %d;X [cm];Y [cm]",
            deid);


        TH2Poly* h =
            new TH2Poly(
                hname,
                htitle,
                xmin,
                xmax,
                ymin,
                ymax);


        // ====================================================
        // Associate histogram with output file
        // ====================================================

        h->SetDirectory(fout);


        // ====================================================
        // Add one polygon/bin for every active pad
        //
        // CURRENT VERSION:
        // rectangular approximation.
        //
        // TODO:
        // replace with TGeo pad polygon.
        // ====================================================

        int nAdded = 0;


        for (const auto& padEntry : pads) {

            const int pid =
                padEntry.first;

            const double charge =
                padEntry.second;


            if (!seg.isValid(pid)) {
                continue;
            }


            // ------------------------------------------------
            // Pad position
            // ------------------------------------------------

            const double x =
                seg.padPositionX(pid);

            const double y =
                seg.padPositionY(pid);


            // ------------------------------------------------
            // Pad dimensions
            // ------------------------------------------------

            const double dx =
                seg.padSizeX(pid) / 2.0;

            const double dy =
                seg.padSizeY(pid) / 2.0;


            const double x1 =
                x - dx;

            const double x2 =
                x + dx;

            const double y1 =
                y - dy;

            const double y2 =
                y + dy;


            // ------------------------------------------------
            // Add rectangular pad
            // ------------------------------------------------

            const int bin =
                h->AddBin(
                    x1,
                    y1,
                    x2,
                    y2);


            if (bin <= 0) {

                std::cerr
                    << "WARNING: could not add pad "
                    << pid
                    << " to DE "
                    << deid
                    << "\n";

                continue;
            }


            // ------------------------------------------------
            // Store ADC
            // ------------------------------------------------

            h->SetBinContent(
                bin,
                charge);


            ++nAdded;
        }


        // ====================================================
        // Write ONLY TH2Poly
        // ====================================================

        fout->cd();

        h->Write();


        std::cout
            << "Valid active pads: "
            << nValidPads
            << "\n";

        std::cout
            << "TH2Poly bins added: "
            << nAdded
            << "\n";


        ++nDE;
    }


    // ========================================================
    // 11. Write output
    // ========================================================

    fout->Write();


    // ========================================================
    // 12. Close files
    // ========================================================

    fout->Close();
    fin->Close();


    delete fout;
    delete fin;


    // ========================================================
    // 13. Final summary
    // ========================================================

    std::cout << "\n";
    std::cout
        << "============================================================\n";

    std::cout
        << "DONE\n";

    std::cout
        << "============================================================\n";

    std::cout
        << "Input file   : "
        << inputFile
        << "\n";

    std::cout
        << "Geometry     : "
        << geometryFile
        << "\n";

    std::cout
        << "Output file  : "
        << outputFile
        << "\n";

    std::cout
        << "Digits       : "
        << nDigits
        << "\n";

    std::cout
        << "DEs          : "
        << nDE
        << "\n";

    std::cout
        << "Invalid pads : "
        << nInvalidPads
        << "\n";

    std::cout << "\n";

    std::cout
        << "Output contains ONLY TH2Poly objects.\n";

    std::cout
        << "No canvas was created or drawn.\n";

    std::cout << "\n";
}