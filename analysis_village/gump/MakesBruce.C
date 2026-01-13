#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>

bool in_array(const std::string &value, const std::vector<std::string> &array)
{
    return std::find(array.begin(), array.end(), value) != array.end();
}

void MakesBruce(const char* fileName = "input.root", const char* output_filename = "output.root") {
    const char* treeNames[1] = {"SelectedEvents"};
    TFile *file = TFile::Open(fileName, "READ");
    TFile *outfile = TFile::Open(output_filename, "recreate");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file '" << fileName << "'" << std::endl;
        return;
    }

    for(int t = 0; t < 1; t++){
        const char* treeName = treeNames[t];
        TTree *tree = nullptr;
        file->GetObject(treeName, tree);
        TObjArray *AllBranches = tree->GetListOfBranches();

        if (!tree) {
            std::cerr << "Error: Could not find TTree '" << treeName << "' in file '" << fileName << "'" << std::endl;
            file->Close();
            delete file;
            return;
        }
        std::cout << "Successfully accessed TTree: " << treeName << std::endl;

        std::map<std::string, std::vector<double>> filler_map;
        std::map<std::string, std::vector<double>> filler_sigmas_map;
        bool is_multisigma = false;
        bool is_multisim = false;
        std::string multisigma_keyword = "multisigma";
        std::string multisim_keyword = "Flux";

        TTree *wgt_multisigma_outtree = new TTree("multisigmaTree", "Systematic weights formatted for PROfit. Using multisigma format.");
        TTree *wgt_multisim_outtree = new TTree("multisimTree", "Systematic weights formatted for PROfit. Using multisim format.");


        for(int b=0; b < AllBranches->GetEntries(); b++){
            TBranch* branch = dynamic_cast<TBranch*>(AllBranches->At(b));
            const char* branchName = branch->GetName();
            std::string branchName_str = branchName;
            std::string branchName_sigma = branchName_str+"_sigma";
            if(branchName_str.find(multisigma_keyword) == std::string::npos && branchName_str.find(multisim_keyword) == std::string::npos){
                continue;
            }

            if(branchName_str.find(multisigma_keyword) != std::string::npos){
                is_multisigma = true;
                std::vector<double> temp(7, 0.0);
                std::vector<double> temp_sigmas(7, 0.0);
                filler_map.insert({branchName_str, temp});
                filler_sigmas_map.insert({branchName_str, temp_sigmas});
                wgt_multisigma_outtree->Branch(branchName, &filler_map[branchName]);
                wgt_multisigma_outtree->Branch(branchName_sigma.c_str(), &filler_sigmas_map[branchName]);
            }

            if(branchName_str.find(multisim_keyword) != std::string::npos){
                is_multisim = true;
                std::vector<double> temp(100, 0.0);
                filler_map.insert({branchName_str, temp});
                wgt_multisim_outtree->Branch(branchName, &filler_map[branchName]);
            }

        }

        Long64_t nEntries = tree->GetEntries();
        double weights_multisigma[7]; 
        double weights_multisim[100];
        for(Long64_t i = 0; i < nEntries; i++){
            for(int b=0; b < AllBranches->GetEntries(); b++){
                TBranch* branch = dynamic_cast<TBranch*>(AllBranches->At(b));
                const char* branchName = branch->GetName();
                std::string branchName_str = branchName;
                std::cout << branchName << std::endl;
                if(branchName_str.find(multisigma_keyword) == std::string::npos && branchName_str.find(multisim_keyword) == std::string::npos){
                    continue;
                }

                 if(branchName_str.find(multisigma_keyword) != std::string::npos){
                    tree->SetBranchAddress(branchName, &weights_multisigma);
                    branch->GetEntry(i);
                    std::vector<double> temp(7, 0.0);
                    std::vector<double> temp_sigmas = {1, -1, 2, -2, 3, -3, 0};
                    for (size_t j = 0; j < 7; j++) {
                        std::cout << "Entry[" << i << "]:" << " Element[" << j << "]: " << weights_multisigma[j] << std::endl;
                        temp.at(j) = weights_multisigma[j];
                    }
                    filler_map[branchName] = temp;
                    filler_sigmas_map[branchName] = temp_sigmas;
                }

                if(branchName_str.find(multisim_keyword) != std::string::npos){
                    tree->SetBranchAddress(branchName, &weights_multisim);
                    branch->GetEntry(i);
                    std::vector<double> temp(100, 0.0);
                    for (size_t j = 0; j < 100; j++) {
                        std::cout << "Entry[" << i << "]:" << " Element[" << j << "]: " << weights_multisim[j] << std::endl;
                        if(isnan(weights_multisim[j])){
                            std::cout << "Found nan!" << std::endl;
                            temp.at(j) = 1.0;
                        }
                        else{
                            temp.at(j) = weights_multisim[j];
                        }
                    }
                    filler_map[branchName] = temp;
                }
                       
                if(i == nEntries - 1){
                    tree->SetBranchStatus(branchName, 0);
                    std::cout << "Done with: " << branchName << std::endl;
                }
            }
           
            wgt_multisigma_outtree->Fill();
            wgt_multisim_outtree->Fill();
        }
            
        tree->CopyTree("");
    }
    std::cout << "--- End of processing ---" << std::endl;
    file->Close();
    delete file;
    outfile->Write();
}
