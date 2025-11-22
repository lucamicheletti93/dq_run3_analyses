/**
 *
 * SIMPLIFIED VERSION OF MACRO
 *
 * This macro generates events with PYTHIA, looks for quarkonia in the events,
 * and stores the information about the events and probes in TTree's that it writes into a root file.
 * The trees can be later used to determine e.g. the multiplicity dependence of J/psi production.
 *
 *
 *
 * General structure:
 *
 *
 *
 *
 *
 * Details:
 *
 * The code was written to be run on a computing farm (like GSI kronos)
 *
 *
 *
 *
 *
 *
 * It takes arguments from the command line:
 *
 * argument 1:  number of events to produce (default: 10000)
 * argument 2:  whether or not to scale down the MB events (0=no, 1=yes) (default: 0)
 * argument 3:  whether to print out more information while running (default: 0)
 * argument 4:  the settings file to use (default: "settings.cmd")
 * argument 5:  the output filename (without the final .root) (default: "out")
 * argument 6:  an initial seed for the random number generator (default:1)
 *
 */



#include <time.h>
#include <vector>
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TTimeStamp.h"
#include "Pythia8/Pythia.h"

//#include "DataTypes.h"


using namespace Pythia8;


/**
 *
 * This is the information stored for each Quarkonium in the tree:
 *  - the pdg code (see http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf )
 *  - the pT and rapidity
 *  - the multiplicities in the towards ("Region1"), transverse ("Region2") and away ("Region3") region
 *    these multiplicities are counted either
 *    - for all etas
 *    - for |eta| < 1
 *    - in the V0A acceptance
 *    - in the V0C acceptance
 *
 */

bool isOnium(int pdg);     // used
bool isPrimaryChargedALICE( unsigned short idx, Pythia &pythia);        //used
bool isLongLived( unsigned int pdg);                                    //used
bool inVZEROacceptance( float eta );                                    //used
bool IsFromJpsi(Particle* part, Pythia& pythia);


/**
 *
 * The main function.
 *
 * - processes the command line arguments
 * - starts parallel threads using openMP
 * - created the output trees
 * - initialized pythia
 * - runs the event loop
 *   - here it calculates the multplicities, looks for the hard probes, and calls helper function
 *     to get the multiplicities in the different regions, and information about the hard probes.
 *
 *
 */
void jpsiHadrons() {
  
  string settings = "settings.cmnd";
  
  int nev = 800000; //7000000
  const char* outfilename = "pythiaJpsiHadron.root";
  
    
    TFile *fout;
    fout = TFile::Open(outfilename,"RECREATE");
    fout->cd();

    // prepare trees
    //EventDQ* event=new EventDQ();
    int eventId, eventType, eventNmpi, eventNtracks, eventMultEta1, eventMultVZERO;
    int oniaPDG, oniaMPDG1, oniaMPDG2;
    float oniaPt, oniaEta, oniaPhi;
    float oniaPx, oniaPy, oniaPz, oniaE;
    float oniaPtLeg1, oniaPtLeg2, oniaEtaLeg1, oniaEtaLeg2, oniaPhiLeg1, oniaPhiLeg2;
    float hadronPt, hadronEta, hadronPhi;
    float hadronPx, hadronPy, hadronPz, hadronE;
    int hadronPDG, hadronPDG1, hadronPDG2;
    bool hadronIsCharged;

    int oniaIdx;
    int oniaD1idx, daughter1PDG;
    int oniaD2idx, daughter2PDG;    

    int batchNo = 1;
    
    TTree* eventTree = new TTree( "eventTree", "event information");
    
    eventTree->Branch("batchNo", &batchNo);
    eventTree->Branch("eventId", &eventId);
    eventTree->Branch("eventType", &eventType);
    eventTree->Branch("eventNmpi", &eventNmpi);
    eventTree->Branch("eventNtracks", &eventNtracks);
    eventTree->Branch("eventMultEta1", &eventMultEta1);
    eventTree->Branch("eventMultVZERO", &eventMultVZERO);
    eventTree->Branch("oniaPDG", &oniaPDG);
    eventTree->Branch("oniaMPDG1", &oniaMPDG1);
    eventTree->Branch("oniaMPDG2", &oniaMPDG2);
    eventTree->Branch("oniaPt", &oniaPt);
    eventTree->Branch("oniaEta", &oniaEta);
    eventTree->Branch("oniaPhi", &oniaPhi);
    eventTree->Branch("oniaPx", &oniaPx);
    eventTree->Branch("oniaPy", &oniaPy);
    eventTree->Branch("oniaPz", &oniaPz);
    eventTree->Branch("oniaE", &oniaE);    
    eventTree->Branch("oniaPtLeg1", &oniaPtLeg1);
    eventTree->Branch("oniaPtLeg2", &oniaPtLeg2);
    eventTree->Branch("oniaEtaLeg1", &oniaEtaLeg1);
    eventTree->Branch("oniaEtaLeg2", &oniaEtaLeg2);
    eventTree->Branch("oniaPhiLeg1", &oniaPhiLeg1);
    eventTree->Branch("oniaPhiLeg2", &oniaPhiLeg2);
    eventTree->Branch("hadronPt", &hadronPt);
    eventTree->Branch("hadronEta", &hadronEta);
    eventTree->Branch("hadronPhi", &hadronPhi);
    eventTree->Branch("hadronPx", &hadronPx);
    eventTree->Branch("hadronPy", &hadronPy);
    eventTree->Branch("hadronPz", &hadronPz);
    eventTree->Branch("hadronE", &hadronE);
    eventTree->Branch("hadronPDG", &hadronPDG);
    eventTree->Branch("hadronPDG1", &hadronPDG1);
    eventTree->Branch("hadronPDG2", &hadronPDG2);
    eventTree->Branch("hadronIsCharged", &hadronIsCharged);
    
    
    // start pythia
    gRandom = new TRandom3();
    Pythia pythia;
    pythia.readFile( settings );

    /*pythia.readString("Random:setSeed = on");
    std::stringstream sstm;
    sstm <<  "Random:seed =" << seed;
    std::string seedString = sstm.str();
    pythia.readString(  seedString  );
    */
    TTimeStamp t;
    int seed = t.GetNanoSec();
    seed = seed % 900000000;    // this is to avoid getting seeds above the pythia threshold

    pythia.readString("Random:setSeed = on");
	std::stringstream sstm;
	sstm <<  "Random:seed =" << seed;
	std::string seedString = sstm.str();
	pythia.readString(seedString);

    pythia.init();

    
    // Start the event loop

    for (int iev = 0; iev < nev; iev++) {

        if (!pythia.next()) continue;
        
        eventId=iev;
        // collision type:
        // single, double, central or non-diffractive, see http://home.thep.lu.se/~torbjorn/pythia82html/QCDProcesses.html
        eventType = pythia.info.code();
        eventNmpi = pythia.info.nMPI();
        eventMultEta1=0; eventNtracks=0; eventMultVZERO=0;
        oniaPDG = 0;
        oniaMPDG1 = 0; oniaMPDG2 = 0; oniaPt = 0.0; oniaEta = 0.0; oniaPhi = 0.0;
        oniaPx = 0; oniaPy = 0; oniaPz = 0; oniaE = 0;
        oniaPtLeg1 = 0.0; oniaEtaLeg1 = 0.0; oniaPhiLeg1 = 0.0;
        oniaPtLeg2 = 0.0; oniaEtaLeg2 = 0.0; oniaPhiLeg2 = 0.0;
        
        oniaIdx = 0;
        
        bool oniumFound=false;
        // run a first loop to check if an onia is present in this event
        for (int iPart = 0;  iPart < pythia.event.size(); iPart++) {
            Particle* part = &pythia.event[iPart];
            if(!oniumFound && 
                abs(part->id())==443 &&
                !IsFromJpsi(part, pythia) &&         
               !(part->id()==pythia.event[part->daughter1()].id()) && 
               !(part->id()==pythia.event[part->daughter2()].id()) && 
               abs(part->y())<0.9) 
            {                                // add kinematic selection on legs, 
                oniumFound = true;
                oniaIdx = iPart;
                oniaPDG = part->id();
                oniaMPDG1 = pythia.event[part->mother1()].id();
                oniaMPDG2 = pythia.event[part->mother2()].id();
                oniaPt = part->pT();
                oniaEta = part->eta();
                oniaPhi = part->phi();
                oniaPx = part->px();
                oniaPy = part->py();
                oniaPz = part->pz();
                oniaE = part->e();
                // at this point we assume the jpsi decay is forced to e+e-, so there are just 2 decay prongs
                oniaPtLeg1 = pythia.event[part->daughter1()].pT();      
                oniaPtLeg2 = pythia.event[part->daughter2()].pT();
                oniaEtaLeg1 = pythia.event[part->daughter1()].eta();
                oniaEtaLeg2 = pythia.event[part->daughter2()].eta();
                oniaPhiLeg1 = pythia.event[part->daughter1()].phi();
                oniaPhiLeg2 = pythia.event[part->daughter2()].phi();
            }
            if(isPrimaryChargedALICE(iPart, pythia)) {
                if(abs(part->eta())<0.9) ++eventNtracks;
                if(abs(part->eta())<1.0) ++eventMultEta1;
                if(inVZEROacceptance(part->eta())) ++eventMultVZERO;
            }
        }
        
        // Skip this event if it does not contain a good jpsi or this event is not selected for the unbiased sample 
        if(!(oniumFound || gRandom->Rndm()<0.0001)) continue;   
        
        cout << "############# Writing event; reason :: ";
        if(oniumFound) cout << "jpsi found!" << endl;
        else cout << "unbiased MB event" << endl;
        cout << "pythia event size :: " << pythia.event.size() << endl;
        cout << "#ALICE primary in eta09/eta1/VZERO :: " << eventNtracks << " / " << eventMultEta1 << " / " << eventMultVZERO << endl;
        
        std::vector<int> daughters;
        std::vector<int> mothers;
        if(oniumFound) {
            std::cout << "jpsi mpdg 1/2 :: " << oniaMPDG1 << " / " << oniaMPDG2 << endl;
            daughters = pythia.event[oniaIdx].daughterList();
            std::cout << "Daughter list for part :" << oniaIdx << "  pdg# " << pythia.event[oniaIdx].id() << endl;
            for(unsigned int i=0; i<daughters.size(); i++) {
                std::cout << "daughter #" << i << "  " << daughters[i] << "  PDG#" << pythia.event[daughters[i]].id() 
                          << "  (pt,eta,phi) : " << pythia.event[daughters[i]].pT() << " / " << pythia.event[daughters[i]].eta() << " / " << pythia.event[daughters[i]].phi() << endl;
            }
            mothers = pythia.event[oniaIdx].motherList();
            std::cout << "Mother list for part :" << oniaIdx << "  pdg# " << pythia.event[oniaIdx].id() << endl;
            for(unsigned int i=0; i<mothers.size(); i++) {
                std::cout << "mother #" << i << "  " << mothers[i] << "  PDG#" << pythia.event[mothers[i]].id() 
                          << "  (pt,eta,phi) : " << pythia.event[mothers[i]].pT() << " / " << pythia.event[mothers[i]].eta() << " / " << pythia.event[mothers[i]].phi() << endl;
            }
        }
        
        cout << "Selecting ALICE primary hadrons: " << endl;
        for(int iPart = 0;  iPart < pythia.event.size(); iPart++) {
            Particle* part = &pythia.event[iPart];
            if(abs(part->y()) > 0.9) continue; // changed from eta to y to avoid break in code!
            if(!isPrimaryChargedALICE(iPart, pythia)) continue;
            if(oniumFound) {
                bool isOniaDaughter = false;
                for(unsigned int i=0; i<daughters.size(); i++)
                    if(iPart==daughters[i])
                        isOniaDaughter = true;
                if(isOniaDaughter) {
                    cout << "track is onia daughter (not added) :: " << iPart << endl;
                    continue;    
                }
            }
            
           // cout << "Adding hadron #" << iPart << " PDG# " << part->id() << endl;
            hadronPt = part->pT();
            hadronEta = part->eta();
            hadronPhi = part->phi();
            
            hadronPx = part->px();
            hadronPy = part->py();
            hadronPz = part->pz();
            hadronE = part->e();
            
            hadronPDG = part->id();
            hadronPDG1 = pythia.event[part->mother1()].id();
            hadronPDG2 = pythia.event[part->mother2()].id();
            
            hadronIsCharged = part->isCharged();
            
            
            eventTree->Fill();
        }
    }

    // here the event loop ends

    eventTree->Write();
    fout->Close();
}

bool inVZEROacceptance( float eta) {
    return (eta > 2.8 && eta < 5.1) || (eta > -3.7 && eta < -1.7);
}

/**
 *
 * Function that finds out, where the quarkonia comes from. In this simplified version it only checks
 * if it is a prompt or non-prompt J/psi (i.e. if one of it's ancestor's was a b-meson or baryon)
 * The function gets the mother of the quarkonium, checks it's id and status code, and decides if it has to
 * go further back to find the initially produced partcile
 * Check http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html for the meaning of the status codes.
 *
 */


// Find out the quark content of a particle from it's PDG identifier
// see http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf


int firstQuark (int pdg) {
    return ( abs(pdg) % 100 - abs(pdg) % 10 ) / 10;
}

int secondQuark (int pdg) {
    return ( abs(pdg) % 1000 - abs(pdg) % 100 ) / 100;
}

int thirdQuark (int pdg) {
    return ( abs(pdg) % 10000 - abs(pdg) % 1000 ) / 1000;
}

bool isCharmonium(int pdg) {
    return ( firstQuark(pdg) == 4 && secondQuark(pdg) == 4 && thirdQuark(pdg) == 0 );
}

bool isBottomonium(int pdg) {
    return ( firstQuark(pdg) == 5 && secondQuark(pdg) == 5 && thirdQuark(pdg) == 0 );
}

bool isOnium(int pdg) {
    return (isCharmonium(pdg) || isBottomonium(pdg) );
}




/**
 *
 * Find out if given partilce is a primary charged particle according to ALICE definition,
 * see https://cds.cern.ch/record/2270008
 *
 */

bool isPrimaryChargedALICE( unsigned short idx, Pythia &pythia ) {

    if( ! pythia.event[idx].isCharged() ) {
        return false;
    }

    int status = pythia.event[idx].statusHepMC();
    if( status==0 || status ==4 || (status>=11 && status<=200) ) return false;
    unsigned int pdg = pythia.event[idx].idAbs();
    if (!isLongLived(pdg)) return false;
    while( (idx = pythia.event[idx].mother1()) ) {
        status = pythia.event[idx].statusHepMC();
        if ( status == 4 || status==13 ) return true;
        pdg = pythia.event[idx].idAbs();
        if ( isLongLived(pdg) ) return false;
    }
    return true;
}


bool isLongLived( unsigned int pdg ) {
    if (pdg > 1000000000) return true;

    switch (pdg) {
    case 13:	//mu⁻
    case 11:	//e⁻
    case 22:	//gamma
    case 211:	//pi⁺
    case 321:	//K⁺
    case 130:	//K⁰_L
    case 310:	//K⁰_S
    case 2212:	//p
    case 2112:	//n
    case 3122:	//lambda
    case 3112:	//sigma⁻
    case 3222:	//sigma⁺
    case 3312:	//xi⁻
    case 3322:	//xi⁰
    case 3334:	//omega⁻
    case 12:	//electron neutrino
    case 14:	//muon neutrino
    case 16:	//tau neutrino
        return true;
    }
    return false;
}

//______________________________________________________________________
bool IsFromJpsi(Particle* part, Pythia& pythia) {
	std::vector<int> mothers = pythia.event[part->index()].motherList();
	for (unsigned int iPart=0; iPart<mothers.size(); iPart++) {
		if (pythia.event[mothers[iPart]].idAbs()==443) return kTRUE;
	}
	return kFALSE;
}
