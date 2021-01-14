from ROOT import gStyle, TCanvas, TH1D, Form, gSystem
gSystem.Load('plottingHelper_C.so')
from ROOT import PlottingHelper as ph



gStyle.SetOptStat(0)

def vec(vv):
    from ROOT import std
    vvv = std.vector("double")()
    for v in vv:
        vvv.push_back(v)
    return vvv


can = TCanvas("can", "")
ph.DividePad(vec([1,1,1,1,1]), vec([1,1,1]))
hh = []
for i in range(5*3):
    can.cd(i+1)
    h = TH1D(""+str(i),";;", 10, -3, 3)
    nEv = 1000 - i*40
    h.FillRandom("gaus", nEv);
    h.Draw("hist e")
    hh.append(h)
    ph.SetFTO(vec([12]), vec([5]), vec([1.2, 2, 0.3, 4]))
    if i == 5*3-1: ph.GetXaxis().SetTitle("x");
    if i == 0: ph.GetYaxis().SetTitle("y");

    ph.GetYaxis().SetRangeUser(0, 300)

    ph.DrawLatexUp( -1,  "n_{Ev} = "+str(nEv));

    #Remove overlaps of both axes
    #RemoveOverlaps(gPad, GetXaxis(), {}, true, true);
    #RemoveOverlaps(gPad, GetYaxis(), {], true, true);

ph.DrawLatexUp(can.GetPad(1), can.GetPad(5), 2, "This is a testing grid in Python");
can.SaveAs("testGridPy.pdf");

