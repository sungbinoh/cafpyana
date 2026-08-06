#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

// This is a root macro for use with the gump/sbruce.py functions for making sbruce files from dataframes.
// Refrain from using this independently, as it doesn't include several preprocessing steps (systematic
// handling, de-duplication, etc.). - Nate, May26

void MakesBruceNew(const char* fileName = "input.root", const char* output_filename = "output.root") {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open input file '" << fileName << "'" << std::endl;
        return;
    }

    TFile *outfile = TFile::Open(output_filename, "RECREATE");
    const char* treeNames[1] = {"SelectedEvents"}; //
    
    for(int t = 0; t < 1; t++) {
        const char* treeName = treeNames[t]; //
        TTree *intree = nullptr;
        file->GetObject(treeName, intree); //

        if (!intree) {
            std::cerr << "Error: Could not find TTree '" << treeName << "' in file '" << fileName << "'" << std::endl; //
            continue;
        }
        std::cout << "Successfully accessed TTree: " << treeName << std::endl; //

        Long64_t nEntries = intree->GetEntries(); //
        TObjArray *allBranches = intree->GetListOfBranches(); //
        int nBranches = allBranches->GetEntries(); //

        // Output trees expected by PROfit
        TTree *wgt_multisigma_outtree = new TTree("multisigmaTree", "Systematic weights formatted for PROfit. Using multisigma format."); //
        TTree *wgt_multisim_outtree = new TTree("multisimTree", "Systematic weights formatted for PROfit. Using multisim format."); //

        // Dynamic vector storage maps
        std::map<std::string, std::vector<double>> vector_storage; //
        std::map<std::string, std::vector<double>> sigma_storage; //

        std::string multisigma_keyword = "multisigma"; //
        std::string multisim_keyword = "Flux"; //

        // 1. Dynamic Branch Booking Loop
        for(int b = 0; b < nBranches; b++) {
            TBranch* branch = dynamic_cast<TBranch*>(allBranches->At(b)); //
            std::string bName = branch->GetName(); //

            if (bName.find(multisigma_keyword) == std::string::npos && bName.find(multisim_keyword) == std::string::npos) { //
                continue;
            }

            TLeaf* leaf = branch->GetLeaf(bName.c_str());
            if (!leaf) continue;

            vector_storage[bName] = std::vector<double>(); //

            if (bName.find(multisigma_keyword) != std::string::npos) { //
                wgt_multisigma_outtree->Branch(bName.c_str(), &vector_storage[bName]); //
                
                std::string sigmaName = bName + "_sigma"; //
                sigma_storage[sigmaName] = std::vector<double>(); //
                wgt_multisigma_outtree->Branch(sigmaName.c_str(), &sigma_storage[sigmaName]); //
            } 
            else if (bName.find(multisim_keyword) != std::string::npos) { //
                wgt_multisim_outtree->Branch(bName.c_str(), &vector_storage[bName]); //
            }
        }

        // 2. Streamlined Event Loop
        for(Long64_t i = 0; i < nEntries; i++) { //
            intree->GetEntry(i); //

            for(int b = 0; b < nBranches; b++) { //
                TBranch* branch = dynamic_cast<TBranch*>(allBranches->At(b)); //
                std::string bName = branch->GetName(); //

                if (vector_storage.find(bName) == vector_storage.end()) continue; //

                TLeaf* leaf = branch->GetLeaf(bName.c_str());
                if (!leaf) continue;

                int arraySize = leaf->GetNdata(); 
                std::vector<double>& current_vec = vector_storage[bName];
                
                // Directly construct elements from the leaf array
                current_vec.assign(arraySize, 0.0);
                for(int j = 0; j < arraySize; j++) {
                    current_vec[j] = leaf->GetValue(j);
                }

                // Assign the required Companion Sigmas based on the discovered array structure
                if (bName.find(multisigma_keyword) != std::string::npos) { //
                    std::string sigmaName = bName + "_sigma"; //
                    std::vector<double>& current_sigma = sigma_storage[sigmaName];
                    
                    if (arraySize == 7) {
                        current_sigma = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0}; //
                    } else if (arraySize == 3) {
                        current_sigma = {0.0, 1.0, -1.0}; //
                    } else if (arraySize == 2) {
                        current_sigma = {0.0, 1.0}; //
                    } else {
                        current_sigma.assign(arraySize, 0.0);
                    }
                }
            }
    	    wgt_multisigma_outtree->Fill(); //
            wgt_multisim_outtree->Fill(); //
        }
        // 3. Metadata Preservation Step
        // Turn off processed weight branches so they aren't duplicated
        for (auto const& [bName, vec] : vector_storage) {
            intree->SetBranchStatus(bName.c_str(), 0); //
        }

        std::cout << "Safely cloning flat kinematic variables into SelectedEvents..." << std::endl;
        
        // Explicitly point ROOT back to the output file so CopyTree writes there directly
        outfile->cd();
        
        // This automatically names the cloned tree "SelectedEvents" matching the input tree
        TTree *clonedTree = intree->CopyTree(""); //
        
        // Re-enable branches on the input tree context in case of further loops
        intree->SetBranchStatus("*", 1);
    } // End of tree loop

    std::cout << "--- End of processing ---" << std::endl; //
    
    // Write out the file, strictly overwriting any historical memory cycles
    outfile->Write(0, TObject::kOverwrite); 
    outfile->Close(); //
    file->Close(); //
    
    delete file; //
    delete outfile;
}
