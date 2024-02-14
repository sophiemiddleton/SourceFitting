// author" S Middleton 2022
#include <fstream>
#include <iostream>
#include <algorithm>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <time.h>

TString filepath = "/Users/user/ski_examples/CaloCalib/10M.root";
using namespace std;
using namespace TMath;

//================== Select the sample =========================//
TTree* get_data_tree(){
    TFile *f =  new TFile(filepath);
    TTree *t = (TTree*)f->Get("CaloExample/Calo");
    return t;
}

int MakeSpec(){
  for(int numVal = 674; numVal < 1348; numVal++){
    std::cout<<" Crystal "<<numVal<<std::endl;

    TString cryNum;
    cryNum = to_string(numVal);

    // create the output file.
    TString ouptName = "mu2e_caloSimu_crySpec" + cryNum + ".root";
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    //create output tree
    Float_t spec, ratio, ntrig, stim, time, tErg;
    TTree *specTree = new TTree("specTree","specTree");
    specTree->Branch("spec", &spec, "spec/F"); //Energy Deposited
    specTree->Branch("ratio", &ratio, "ratio/F"); //Ratio of Energy Deposited
    specTree->Branch("ntrig", &ntrig,"ntrig/F");//Number of Hits on Target Crystal
    specTree->Branch("stim", &stim, "stim/F");//Time Difference from Same Crystal
    specTree->Branch("time", &time, "time/F");//Maxium Time Difference for All Hits
    specTree->Branch("tErg", &tErg, "tErg/F");//Energy Spectrum with Special Time Difference

    // input tree
    TTree *ergTree = get_data_tree();

    float cryEdep[100];
    int nCry;
    int cryId[100];
    float cryTime[100];
    float cryPosX[100], cryPosY[100], cryPosZ[100];

    ergTree -> SetBranchAddress("nCry", &nCry);
    ergTree -> SetBranchAddress("cryId", &cryId);
    ergTree -> SetBranchAddress("cryEdep", &cryEdep);
    ergTree -> SetBranchAddress("cryTime", &cryTime);
    ergTree -> SetBranchAddress("cryPosX", &cryPosX);
    ergTree -> SetBranchAddress("cryPosY", &cryPosY);
    ergTree -> SetBranchAddress("cryPosZ", &cryPosZ);

    unsigned int nEvt = (int)ergTree -> GetEntries();
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      ergTree -> GetEntry(iEvt);
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(cryEdep[icry]);
        }
      }
      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float sameCryDeltaT = sameCryTime.back() - sameCryTime.front();
      stim = sameCryDeltaT;

      float ergTarget = 0.0;
      float ergRest   = 0.0;
      int nTarget = 0;
      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          ergTarget = cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
          nTarget += 1;
        }
        else
        {
          ergRest += cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      int numHit = sameCryEdep.size();
      if(difTime > 20)
      {
        for(int iHit = 0; iHit < numHit; iHit++)
        {
          tErg = sameCryEdep[iHit];
        }
      }
      if((ergTarget / (ergTarget + ergRest)) >= 0.8 && difTime < 4)
      {
        spec =ergTarget;
      }
      ratio= (ergTarget / (ergTarget + ergRest));
      ntrig = (nTarget);
      time  = (difTime);
      specTree->Fill();
    }

  ouptFile -> Write();
  ouptFile -> Close();

  }
    return 0;
}
