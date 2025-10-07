#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>

// This macro reads a branch containing a vector of doubles from a TTree
// and prints the elements of the vector for each entry.

void MakesBruce(const char* fileName = "input.root", const char* output_filename = "output.root") {
    // --- Configuration ---
    // Replace "your_file.root" with the actual name of your ROOT file.
    const char* treeName = "SelectedEvents";

    // --- Open the ROOT file ---
    // The "READ" option opens the file in read-only mode.
    TFile *file = TFile::Open(fileName, "READ");

    // Check if the file was opened successfully.
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file '" << fileName << "'" << std::endl;
        return;
    }
    std::cout << "Successfully opened file: " << fileName << std::endl;

    // --- Get the TTree ---
    TTree *tree = nullptr;
    file->GetObject(treeName, tree);
    TObjArray *AllBranches = tree->GetListOfBranches();

    // Check if the TTree was retrieved successfully.
    if (!tree) {
        std::cerr << "Error: Could not find TTree '" << treeName << "' in file '" << fileName << "'" << std::endl;
        file->Close();
        delete file;
        return;
    }
    std::cout << "Successfully accessed TTree: " << treeName << std::endl;

    TFile *outfile = TFile::Open(output_filename, "recreate");
    TTree *wgt_outtree = new TTree("multisigmaTree", "Systematic weights formatted for PROfit");

    tree->SetBranchStatus("*", 0);
    for(int b=0; b < AllBranches->GetEntries(); b++){
        TBranch* branch = dynamic_cast<TBranch*>(AllBranches->At(b));
        const char* branchName = branch->GetName();
        std::cout << b << ", " << branchName << std::endl;

        std::string branchName_str = branchName;
        std::string keyword = "multisigma";
        if(branchName_str.find(keyword) == std::string::npos){
            std::cout << "Not treating as weight because it doesn't contain the keyword: " << keyword << std::endl;
            tree->SetBranchStatus(branchName, 1);
            continue;
        }

        // --- Set up the branch reading ---
        // We assume the branch contains std::vector<double>.
        // A pointer to a vector of doubles is created to hold the data.
        double weights[7];

        // Link the branch in the TTree to our C++ vector pointer.
        // ROOT will handle the memory allocation for the vector.
        tree->SetBranchAddress(branchName, &weights);

        // --- Loop over all entries in the TTree ---
        Long64_t nEntries = tree->GetEntries();
        std::cout << "Processing " << nEntries << " entries..." << std::endl;
        std::vector<double> event_wgts;
        wgt_outtree->Branch(branchName, &event_wgts);
        for (Long64_t i = 0; i < nEntries; ++i) {
            // Load the data for the i-th entry. This fills our 'weights' vector.
            tree->GetEntry(i);
            event_wgts.clear();

            // Check if the weights vector is valid and not empty for this entry.
            // Loop through the elements of the vector and print them.
            for (size_t j = 0; j < 7; ++j) {
                std::cout << "  Element[" << j << "]: " << weights[j] << std::endl;
                event_wgts.push_back(weights[j]);
            }
            wgt_outtree->Fill();
        }

    }
    tree->CopyTree("");
    // --- Clean up ---
    std::cout << "--- End of processing ---" << std::endl;
    file->Close();
    delete file; // Good practice to free the memory.
    outfile->Write();
}
