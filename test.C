#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TGaxis.h"

#include "plottingHelper.h"
#include "RemoveOverlaps.h"

R__LOAD_LIBRARY(plottingHelper_C.so)
using namespace PlottingHelper;


void test()
{
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "");
    DividePad({1,1,1,1,1}, {1,1,1});
    for(int i = 0; i < 5*3; ++i) {
        can->cd(i+1);
        
        TH1D *h = new TH1D(Form("%d",rand()),";;", 10, -3, 3);
        int nEv = 1000 - i*40;
        h->FillRandom("gaus", nEv);
        h->Draw("hist e");
        SetFTO({12}, {5}, {1.2, 2, 0.3, 4});
        if(i == 5*3-1) GetXaxis()->SetTitle("x");
        if(i == 0) GetYaxis()->SetTitle("y");
        GetYaxis()->SetRangeUser(0, 300);
        DrawLatexUp( -1,  Form("n_{Ev} = %d", nEv));
        //Remove overlaps of both axes
        RemoveOverlaps(gPad, GetXaxis(), {"-1", "2", "3"}, true, true);
        RemoveOverlaps(gPad, GetYaxis(), {"150", "200"}, true, true);
    }

    DrawLatexUp(can->GetPad(1), can->GetPad(5), 2, "This is a testing grid");
    can->SaveAs("testGrid.pdf");

}

int main()
{
    test();
    return 0;
}   
