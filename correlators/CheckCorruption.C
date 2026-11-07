// Scan every z-auau_run_*.root file individually (TChain silently skips bad
// files, so per-file checks are required) and report which ones are corrupt.
// A file is flagged bad if: TFile::Open fails / IsZombie, the "teposevent"
// tree is missing, GetEntries()==0, or reading the last entry's branches
// throws/fails (catches truncated mid-write files with a valid header but
// damaged tail baskets).
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TError.h>
#include <fstream>
#include <vector>
#include <string>

void CheckCorruption() {
    gErrorIgnoreLevel = kFatal; // suppress ROOT's own recovery warnings; we do our own check

    TSystemDirectory dir(".", ".");
    TList *files = dir.GetListOfFiles();
    std::vector<std::string> names;
    if (files) {
        TSystemFile *f;
        TIter next(files);
        while ((f = (TSystemFile*)next())) {
            std::string name = f->GetName();
            if (!f->IsDirectory() && name.rfind("z-auau_run_", 0) == 0 &&
                name.size() > 5 && name.substr(name.size()-5) == ".root") {
                names.push_back(name);
            }
        }
    }
    std::sort(names.begin(), names.end());

    std::ofstream bad("bad_files.txt");
    std::ofstream good("good_files.txt");
    int nbad = 0, ngood = 0;
    long long totalEntries = 0;

    for (size_t i = 0; i < names.size(); ++i) {
        const std::string &fn = names[i];
        TFile *tf = TFile::Open(fn.c_str(), "READ");
        bool isBad = false;
        std::string reason;

        if (!tf || tf->IsZombie()) {
            isBad = true;
            reason = "open_failed_or_zombie";
        } else {
            TTree *tree = (TTree*)tf->Get("teposevent");
            if (!tree) {
                isBad = true;
                reason = "no_teposevent_tree";
            } else {
                Long64_t n = tree->GetEntries();
                if (n <= 0) {
                    isBad = true;
                    reason = "zero_entries";
                } else {
                    // try reading last entry to catch truncated tail baskets
                    tree->SetBranchStatus("*", 0);
                    tree->SetBranchStatus("np", 1);
                    tree->SetBranchStatus("bim", 1);
                    Int_t nb = tree->GetEntry(n - 1);
                    if (nb <= 0) {
                        isBad = true;
                        reason = "last_entry_read_failed";
                    } else {
                        totalEntries += n;
                    }
                }
            }
        }

        if (isBad) {
            printf("BAD  %-30s %s\n", fn.c_str(), reason.c_str());
            bad << fn << "\n";
            nbad++;
        } else {
            good << fn << "\n";
            ngood++;
        }

        if (tf) { tf->Close(); delete tf; }
    }

    printf("\n=== Summary: %d good, %d bad, %zu total ===\n", ngood, nbad, names.size());
    printf("Total events in good files: %lld\n", totalEntries);
}
