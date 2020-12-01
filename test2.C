#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TGaxis.h"

#include "plottingHelper.h"
#include "RemoveOverlaps.h"

R__LOAD_LIBRARY(plottingHelper_C.so)
using namespace PlottingHelper;


void test2()
{
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "", 1200, 800);
    DivideTransparent({1,0,1,0,1,0,1,0,1}, {1,0,1,0,1});
    for(int i = 0; i < 5*3; ++i) {
        can->cd(i+1);
        
        TH1D *h = new TH1D(Form("%d",rand()),";;", 100, 30, 50000);
        int nEv = 10*(1000 - i*40);
        for(int j = 0; j < nEv; ++j) {
            double v = exp(pow(gRandom->Uniform()*10,2));
            h->Fill(v);
        }
        h->Draw("hist e");
        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();



        gPad->SetLogx();
        SetFTO({12}, {5}, {1.4, 2, 0.5, 4});
        if(i == 5*3-1) GetXaxis()->SetTitle("x");
        if(i == 0) GetYaxis()->SetTitle("y");
        GetYaxis()->SetRangeUser(0, 300);
        DrawLatexUp( -1,  Form("n_{Ev} = %d", nEv));
        //Remove overlaps of both axes

        if(i%5 != 0 ) GetYaxis()->SetLabelSize(0.001);
        if(i < 10 ) GetXaxis()->SetLabelSize(0.001);

        RemoveOverlaps(gPad, GetXaxis(), true, true);
        RemoveOverlaps(gPad, GetYaxis(), true, true);
    }

    DrawLatexUp(can->GetPad(1), can->GetPad(5), 2, "This is a testing grid");
    can->SaveAs("testGrid.pdf");

}

int main()
{
    test2();
    return 0;
}   
