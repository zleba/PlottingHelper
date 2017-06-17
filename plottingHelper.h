#ifndef __plottingHelper__
#define __plottingHelper__

#include "TPad.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>
#include "TString.h"
#include "TFrame.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TMarker.h"
#include "THStack.h"
#include "TMultiGraph.h"

#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <string>

/// The namespace of whole Plotting Helper utility
/// 
///
namespace plottingHelper {

using namespace std;

#define SF TString::Format 



/// Contains coordinates object rectangle envelope
struct Rectangle {
	double fX, fY, fWidth, fHeight;
};


struct Point {
	Point(): x(0), y(0) {}
	Point(double _x, double _y): x(_x), y(_y) {}
	double x, y;
};



/// Enum with positions which are used for legend placing.
///
/// The numbers corresponds to numeric keyboard on Nokia cell-phones
/// It means 1 - top left, 3 top right
/// 7 - bottom left, 9 bottom right
/// the side positions 2,4,6,8 denotes that legend can be placed on arbitrary place along the side
/// The kPos5 represents arbitrary position
enum {
    kPos1 = BIT(1), ///< Left top corner
    kPos2 = BIT(2), ///< Top side
    kPos3 = BIT(3), ///< Right top corner
    kPos4 = BIT(4), ///< Left side
    kPos5 = BIT(5), ///< Arbitrary position
    kPos6 = BIT(6), ///< Right side
    kPos7 = BIT(7), ///< Left bottom corner
    kPos8 = BIT(8), ///< Bottom side
    kPos9 = BIT(9)  ///< Right bottom corner
};



/// Iterator which goes over all positions specified in the constructor
///
/// Note that bitwise "or" operator | can be used to select more complex layout
///
struct PosIterator {
    PosIterator(int _nSteps, unsigned Pos) {
        pos = Pos;
        nSteps = _nSteps;
        last = nSteps-1;
        iSave = -1; jSave = nSteps;
        //Iterate();
    }

    bool Iterate() {
        int last = nSteps-1;
        int i=0, j=0;
        //cout << iSave <<" "<< jSave <<" "<< nSteps<< endl;
        for(i = iSave; i < nSteps; ++i)
        for(j = (jSave+1)*(i==iSave); j < nSteps; ++j) {
            if(pos & kPos5)
                goto after;
            if((pos & kPos3)  && i==0 && j==last)
                goto after;
            if((pos & kPos1)  && i==0 && j==0)
                goto after;
            if((pos & kPos7)  && i==last && j==0)
                goto after;
            if((pos & kPos9)  && i==last && j==last)
                goto after;
            if((pos & kPos2)  && i==0)
                goto after;
            if((pos & kPos8)  && i==last)
                goto after;
            if((pos & kPos4)  && j==0)
                goto after;
            if((pos & kPos6)  && j==last)
                goto after;

        }
        after:
        iSave = i;
        jSave = j;
        //cout << "Radek " << i << " "<< j << endl;
        if(i > last)
            return false;
        else
            return true;
    }

    int nSteps, last;
    int iSave, jSave;
    unsigned pos;
};

/// @name Getters
/// Functions returning elements of the active frame
///@{ 

/// The method returns current frame on the active pad.
///
/// Useful for example to set maximum of the y-axis, GetFrame()->SetMaximum(5)
inline TH1 *GetFrame()
{
   // get first histogram in the list of primitives

	TH1 *hobj = nullptr;

   TIter next(gPad->GetListOfPrimitives());
   TObject *obj;
   while ((obj = next())) {
      if(obj->InheritsFrom(TH1::Class())) {
         hobj = (TH1*)obj;
         //hobj->DrawCopy("sameaxis");
			break;
      }
      if(obj->InheritsFrom(TMultiGraph::Class())) {
         TMultiGraph *mg = (TMultiGraph*)obj;
			hobj = mg ? mg->GetHistogram() : nullptr;
			//if (hobj) hobj->DrawCopy("sameaxis");
			break;
      }
      if(obj->InheritsFrom(TGraph::Class())) {
         TGraph *g = (TGraph*)obj;
			hobj = g ? g->GetHistogram() : nullptr;
         //hobj->DrawCopy("sameaxis");
			break;
      }
      if(obj->InheritsFrom(THStack::Class())) {
         THStack *hs = (THStack*)obj;
			hobj = hs ? hs->GetHistogram() : nullptr;
			//hobj->DrawCopy("sameaxis");
			break;
      }
   }
	return hobj;
}

/// Return x-axis of the active frame
inline TAxis *GetXaxis() { return GetFrame()->GetXaxis(); }  

/// Return y-axis of the active frame
inline TAxis *GetYaxis() { return GetFrame()->GetYaxis(); }  

/// Return z-axis of the active frame
inline TAxis *GetZaxis() { return GetFrame()->GetZaxis(); }  

///@} 



/// @name Units helpers
/// Functions to convert units and obtain axis fractions
///@{ 

/// Transform font size from px to the relative units
///
/// The returned size is dependent on the currently active pad - the gPad pointer
/// @param px The font size in px
/// @return The font size in relative units
inline double PxFontToRel(double px)
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();

	if(pxW < pxH) return px/pxW;
	         else return px/pxH;

    //double pad_width  = gPad->XtoPixel(gPad->GetX2());
    //double pad_height = gPad->YtoPixel(gPad->GetY1()); 
	//if(pxW < pxH) return px/pad_width / corX;
	         //else return px/pad_height / corY;
}


/// Transform font size from relative units to pixels
///
/// The returned size is dependent on the currently active pad - the gPad pointer
/// @param rel The font size in relative units
/// @return The font size in pixels
inline double RelFontToPx(double rel)
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();
	if(pxW < pxH) return rel*pxW;
	         else return rel*pxH;
}


/// Fraction of the pad width covered by x-axis
inline double GetAxisFractionX()
{
    //return (gPad->GetUxmax() - gPad->GetUxmin())/(gPad->GetX2() - gPad->GetX1());
    return (1. - gPad->GetLeftMargin() - gPad->GetRightMargin());
}

/// Fraction of the pad height covered by y-axis
inline double GetAxisFractionY()
{
    //return (gPad->GetUymax() - gPad->GetUymin())/(gPad->GetY2() - gPad->GetY1());
    return (1. - gPad->GetBottomMargin() - gPad->GetTopMargin());
}

/// Convert tick length in px of x-axis to the relative units
inline double TickAbsToRelX(double tick)
{
    double nom = GetAxisFractionX() * gPad->GetAbsHNDC() * gPad->GetWh();
	return tick/nom;
}


/// Convert tick length in px of y-axis to the relative units
inline double TickAbsToRelY(double tick)
{
    double nom =  GetAxisFractionY() * gPad->GetAbsWNDC() * gPad->GetWw();
	return tick/nom;
}

///@}

/// @name Decorators
/// Functions to simply define text sizes, offsets and axis ticks
///@{ 

/// Set text size of x- and y-axis titles and labels to the frame title
///
/// Note that the Canvas width and height is specified in TCanvas constructor.
/// The font size is defined as the height of envelope between "gh" characters
/// The height of numbers is for example 3/4 of the font size.
inline void SetFonts(double pxX, double pxY = -1, double pxT = -1)
{
    TH1 *frame = GetFrame();
    if(pxY < 0) pxY = pxX;
    if(pxT < 0) pxT = pxX;

    double relSizeX = PxFontToRel(pxX);
    double relSizeY = PxFontToRel(pxY);
    double relSizeT = PxFontToRel(pxT);

	frame->GetXaxis()->SetTitleSize(relSizeX);
	frame->GetXaxis()->SetLabelSize(relSizeX);

	frame->GetYaxis()->SetTitleSize(relSizeY);
	frame->GetYaxis()->SetLabelSize(relSizeY);

    frame->SetTitleSize(relSizeT);
}

///
/// Set ticks sizes of x and y-axis in px
///
inline void SetTicks(double tickX, double tickY = -1)
{
    //gPad->Update();
    if(tickY < 0) tickY = tickX;
    TH1 *frame = GetFrame();
	frame->GetXaxis()->SetTickSize(TickAbsToRelX(tickX) );
	frame->GetYaxis()->SetTickSize(TickAbsToRelY(tickY) );
}

/// Set offset of the x-axis labels.
///
/// The zero offset correspond to labels "standing" on top of the x-axis.
/// The offset is measured in units corresponds to height of numbers like "012..."
/// Therefore offset 1 correspond to labels touching x-axis by from the bottom.
inline void SetLabelOffsetX(double off)
{
    TH1 *frame = GetFrame();
    double labSize = frame->GetXaxis()->GetLabelSize();
    double off0 = 0;
    if(gPad->GetLogx()) {
        off0 = -1.15;
        if(frame->GetXaxis()->GetNoExponent()) off0 += 0.33;
    }
    else {
        off0 = -0.8;
    }

    double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
    double fact = RelFontToPx(labSize)/pxH;
    off *= fact *3/4.; //conversion to "1354" height


    off += off0 * labSize;
    frame->GetXaxis()->SetLabelOffset(off);
}




/// Set offset of the x-axis title.
///
/// The zero offset corresponds to title which is intersect by the x-axi in its middle.
/// The offset is measured in units corresponds to height of numbers like "012..."
/// Therefore offset 0.5 corresponds to numeric title touching x-axis by from the bottom.
/// Note that the height of numbers is 3/4 of the font size (height) which is defined by the
/// gh vertical envelope. It means that if title contains h, offset 4/3 * 0.5 would result of "h" touching axis
/// The title text size must be set before offset!
///
inline void SetTitleOffsetX(double off)
{
    TH1 *frame = GetFrame();
    TAxis *ax = frame->GetXaxis();
    off *= 3/4.;
    ax->SetTitleOffset(off* RelFontToPx(ax->GetTitleSize())/1.6 / ( gPad->GetWh()* ax->GetTitleSize()*gPad->GetAbsHNDC() ) );
}

/// Set offset of the y-axis labels.
///
/// The zero offset corresponds to labels touching the y-axis from the left
/// The offset is measured in units corresponds to height of numbers like "012..."
/// I.e. the width of 0 if lying.
/// The label text size must be set before offset!
///
inline void SetLabelOffsetY(double off)
{
    TH1 *frame = GetFrame();
    TAxis *ax = frame->GetYaxis();
    double labSize = ax->GetLabelSize();
    double pxW =  gPad->GetWw() * gPad->GetAbsWNDC() ;
    double fact = RelFontToPx(labSize)/pxW;

    double off0 = 0;
    if(gPad->GetLogy())
        off0 = -0.08;


    off *= 3/4.;
    ax->SetLabelOffset(off *  fact + off0*labSize );
}

/// Set offset of the y-axis title.
///
/// The zero offset corresponds to title which is intersect by the y-axis in its middle.
/// The offset is measured in units corresponds to height of numbers like "012..."
/// Note that the height of numbers is 3/4 of the font size (height) which is defined by the
/// gh vertical envelope.
/// The title text size must be set before offset!
///
inline void SetTitleOffsetY(double off)
{
    TH1 *frame = GetFrame();
    TAxis *ax = frame->GetYaxis();
    off *= 3/4.;//to vertical size of "0"
    ax->SetTitleOffset( off*  RelFontToPx(ax->GetTitleSize())/1.6 / ( gPad->GetWw()*ax->GetTitleSize()* gPad->GetAbsWNDC() ) );
}

/// Set of the fonts and ticks sizes 
///
/// Set size of x-,y-axis labels and titles to px (in pixels).
/// In addition the x and y tick lengths can be specified
inline void SetFontsTicks(double px, double tickX, double tickY = -1.)
{
    SetFonts(px);
    SetTicks(tickX, tickY);
}

/// Set offsets of x,y labels and titles
///
/// Note that root is using different alignment of labels and titles
/// We don't try to correct for that as we assume that in general
/// the actual form of the title can be specified later
/// The title's and label's text sizes must be set before offsets!
inline void SetOffsets(double lX, double tX, double lY, double tY)
{
    SetLabelOffsetX(lX);
    SetTitleOffsetX(tX);
    SetLabelOffsetY(lY);
    SetTitleOffsetY(tY);
}

/// Set Fonts, Ticks and Offsets
///
/// Set up fonts sizes, ticks length and titles/labels offsets
/// All parameters are given as vectors, call as SetFTO({14}, {5}, {1.3, 2.3, 0.3, 2.3});
inline void SetFTO(vector<double> fonts, vector<double> ticks, vector<double> offsets)
{
    fonts.resize(3, -1);
    ticks.resize(2, -1);
    offsets.resize(4, 0);

    SetFonts(fonts[0], fonts[1], fonts[2]);
    SetTicks(ticks[0], ticks[1]);
    SetOffsets(offsets[0], offsets[1], offsets[2], offsets[3]);
}

///@}



/// @name General titles
/// Functions to simplify drawing of pad captions outside of the frames
/// Especially useful in case of describing complex grid of frames.
/// For captions inside of frame consider using of automatic legend.
/// Currently under development.
///@{ 


/// Draw text in top of the current pad.
///
/// Text size is given by frame title
/// The text alignment l,c,r can be specified
///
inline void DrawTextUp(TString text, TString pos = "c")
{
    double t = gPad->GetTopMargin();
    double l = gPad->GetLeftMargin();
    double r = gPad->GetRightMargin();

    double fSize = GetFrame()->GetTitleSize();

    double off = 0.7;
    double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
    double fact = RelFontToPx(fSize)/pxH;
    off *= fact *3/4.; //conversion to "1354" height

    TLatex *tex = new TLatex;
    tex->SetTextSize(fSize);
    if(pos == "c") {
        tex->SetTextAlign(21);
        tex->DrawTextNDC((l+1-r)/2, 1-t+off, text);
    }
    else if(pos == "l") {
        tex->SetTextAlign(11);
        tex->DrawTextNDC(l, 1-t+off, text);
    }
    else if(pos == "r") {
        tex->SetTextAlign(31);
        tex->DrawTextNDC(1-r, 1-t+off, text);
    }
}

inline TLatex *CreateText(TString text, double fSize, int align)
{
    TVirtualPad *padOrg = gPad;
    gPad->GetCanvas()->cd();
    TLatex *lat = new TLatex(0, 0, text);
    lat->SetNDC();
    lat->SetTextSize(PxFontToRel(fSize));
    lat->SetTextAlign(align);
    padOrg->cd();
    return lat;
}


inline void DrawText(TVirtualPad *pad1, TVirtualPad *pad2, TLatex *lat, unsigned pos, double dist)
{
    TVirtualPad *padOrg = gPad;
    assert(pad1->GetCanvas() == pad2->GetCanvas());
    TCanvas *can = pad1->GetCanvas();
    can->cd();

    double l[2], r[2], b[2], t[2];

    TVirtualPad *pads[] = {pad1, pad2};

    bool isHorisontal = true;
    if(pad1->GetAbsYlowNDC() == pad2->GetAbsYlowNDC())
        isHorisontal = true;
    else if(pad1->GetAbsXlowNDC() == pad2->GetAbsXlowNDC())
        isHorisontal = false;
    else
        return;

    for(int i = 0; i < 2; ++i) {
        l[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*pads[i]->GetLeftMargin();
        r[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*(1-pads[i]->GetRightMargin());
        b[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*pads[i]->GetBottomMargin();
        t[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*(1-pads[i]->GetTopMargin());
    }
    double left  = min(l[0], l[1]);
    double right = max(r[0], r[1]);
    cout << "LeftRight " << left <<" "<< right << endl;
    double bottom= min(b[0], b[1]);
    double top   = max(t[0], t[1]);

    double x, y;
    if(isHorisontal) {
        if     (pos == kPos1 || pos == kPos7) x = left;
        else if(pos == kPos2 || pos == kPos8) x = (left + right) / 2;
        else if(pos == kPos3 || pos == kPos9) x = right;
        else return;

        double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
        double fact = RelFontToPx(lat->GetTextSize())/pxH;
        dist *= fact * 3/4.; //conversion to "1354" height


        if(pos == kPos7 || pos == kPos8 || pos == kPos9)
            y = pad1->GetAbsYlowNDC() - dist;
        else if(pos == kPos1 || pos == kPos2 || pos == kPos3)
            y = pad1->GetAbsYlowNDC() + (1-pad1->GetTopMargin())*pad1->GetAbsHNDC() + dist;
        else return;
        cout <<"yDist "<< pad1->GetAbsYlowNDC() <<" "<< pad1->GetAbsHNDC()<<" " << dist <<endl;
    }
    else { //verticaly located pads
        cout << "Is vertical" << endl;
        if     (pos == kPos7 || pos == kPos9) y = bottom;
        else if(pos == kPos4 || pos == kPos6) y = (top + bottom) / 2;
        else if(pos == kPos1 || pos == kPos3) y = top;
        else return;

        double pxW =  gPad->GetWw() * gPad->GetAbsWNDC();
        double fact = RelFontToPx(lat->GetTextSize())/pxW;
        dist *= fact * 3/4.;

        cout <<"xDist "<< pad1->GetAbsXlowNDC() <<" "<< pad1->GetAbsWNDC()<<" " << dist <<endl;
        if(pos == kPos1 || pos == kPos4 || pos == kPos7)
            x = pad1->GetAbsXlowNDC() - dist;
        else if(pos == kPos3 || pos == kPos6 || pos == kPos9)
            x = pad1->GetAbsXlowNDC() + pad1->GetAbsWNDC() + dist;
        else return;
    }

    if(pos == kPos1 || pos == kPos2 || pos == kPos3) {
        //double l1 = pad1->GetAbsXlowNDC() + pad1->GetAbsWNDC()*pad1->GetLeftMargin();
        //double l2 = pad2->GetAbsXlowNDC() + pad2->GetAbsWNDC()*pad2->GetLeftMargin();
    }

    cout << "Plotting " <<(void*)gPad<<" "<< x <<" "<< y<<" "<< lat->GetTitle() << endl;
    lat->DrawTextNDC(x, y, lat->GetTitle());

    gPad->Update();
    padOrg->cd();
}


///@}







///
/// Struct specifying borders of some object which can be approximated by rectangles
/// We use such approach to define position of histograms
///
struct Borders{
    vector<Rectangle> recs;
    void GetHistBorders(TH1 *h);
    void FromAbs2px(double scaleUp, double scaleDn);
    void FromNDC2px();
    static double Distance2(const Rectangle &r1, const Rectangle &r2);
};

//Don't search for lower distance if minSkip reach
static double MinDistance2(double minSkip, const Borders &br1, const Borders &br2);

//Don't search for lower distance if minSkip reach
static double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h);
static double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip);





static Point Px2NDC(Point p);

struct PLACER {

    int nLeg;
    vector<unsigned> dims;
    vector<vector<pair<int,int>>> layOuts;
    vector<vector<double>> distToHists;

    vector<vector<Borders>> legBorders;
    vector<Borders> borders;

    vector<double> SizesX, SizesY;

    int nSteps = 10;
    double minSepar;

    //Initalize layOuts, dims, legBorders, borders
    void init(vector<TLegend *> legs, vector<double> &sizesX, vector<double> &sizesY) {
        nLeg = legs.size();
        for(auto & leg : legs) {
            unsigned layOut = atoi(leg->GetName());

            vector<pair<int,int>> pos;

            PosIterator iter(nSteps, layOut);
            while(iter.Iterate()) {
               pos.push_back(make_pair(iter.iSave, iter.jSave));
            }
            layOuts.push_back(pos);
        }

        dims.resize(layOuts.size());
        for(unsigned i = 0; i < dims.size(); ++i)
            dims[i] = layOuts[i].size();
        SizesX = sizesX;
        SizesY = sizesY;
        GetLegendsPositions();

    }

    double iterate(vector<int> &bestLayout) {
        vector<int> idxs(dims.size());
        vector<double> distsNow(dims.size());
        vector<double> distsBest(dims.size());

        double maxDist = 0;

        while (1) {
            // Print
            double dist = analyze(idxs, distsNow);

            if(dist >= maxDist) {
                bool state = true;
                for(int i = 0; i < nLeg; ++i) {
                    state = state && distsNow[i] >= distsBest[i];
                }
                if(dist > maxDist || (dist == maxDist && state)) {
                    bestLayout = idxs;
                    distsBest = distsNow;
                }
                maxDist = dist;
            }

            // Update
            unsigned j;
            for(j = 0; j < dims.size(); ++j) {
                idxs[j]++;
                if (idxs[j] < static_cast<int>(dims[j])) break;
                idxs[j] = 0;
            }
            if (j == dims.size()) break;
        }
        return maxDist;
    }


    void GetLegendsPositions() {

        LoadHistoBorders(borders);

        legBorders.resize(nLeg);

        double fontSize = RelFontToPx(GetYaxis()->GetTitleSize());
        minSepar = 0.5*fontSize; //inPixels

        TH1 *frameNow = GetFrame();
        double ytickLowNDC = frameNow->GetYaxis()->GetTickLength() * GetAxisFractionY();
        double ytickHighNDC = ytickLowNDC * gPad->GetTicky();
        double xtickLowNDC = frameNow->GetXaxis()->GetTickLength() * GetAxisFractionX();
        double xtickHighNDC = xtickLowNDC * gPad->GetTickx();

        double fSizeNDCy = minSepar / (gPad->GetWh() * gPad->GetAbsHNDC());
        double fSizeNDCx = minSepar / (gPad->GetWw() * gPad->GetAbsWNDC());


        for(int k = 0; k < nLeg; ++k) {

            double SizeX = SizesX[k];
            double SizeY = SizesY[k];

            double xMin =  gPad->GetLeftMargin() + ytickLowNDC + fSizeNDCx;
            double xMax =  1 - gPad->GetRightMargin() - ytickHighNDC - SizeX - fSizeNDCx;
            double yMax =  1 - gPad->GetTopMargin() - xtickHighNDC - SizeY - fSizeNDCy;
            double yMin =  gPad->GetBottomMargin() + xtickLowNDC  + fSizeNDCy;


            legBorders[k].resize(dims[k]);

            for(unsigned l = 0; l < dims[k]; ++l) {
                int i = layOuts[k][l].first;
                int j = layOuts[k][l].second;
                
                double xx = xMin + (xMax-xMin)/(nSteps-1) * j;
                double yy = yMax - (yMax-yMin)/(nSteps-1) * i;


                //Legend borders
                legBorders[k][l].recs.push_back({xx, yy, SizeX, SizeY});
                legBorders[k][l].FromNDC2px();

            }
        }
    }


    void GetDistancesPx(double scaleUp, double scaleDn) {
        auto bordersPx = borders;
        for(auto &b : bordersPx)
            b.FromAbs2px(scaleUp, scaleDn);

        distToHists.resize(nLeg);

        for(int i = 0; i < nLeg; ++i) {
            distToHists[i].resize(dims[i]);
            for(unsigned l = 0; l < dims[i]; ++l) 
                distToHists[i][l] = MinDistanceSingle(bordersPx, legBorders[i][l], 0.99*minSepar);
        }
    }


    pair<double,double> GetSolution(vector<double> &xx, vector<double> &yy,
                                                         int nScaleSteps=10) {
        
        double scaleUp=1., scaleDn=1.;

        vector<int> bestLayout;

        for(int sum = 0; sum < nScaleSteps; ++sum)
        for(int iUp = sum; iUp >= 0; --iUp) {
                int iDn = sum - iUp;
                if(!gPad->GetLogy() && iDn != 0) continue;
                scaleUp = pow(4, iUp/10.);
                scaleDn = pow(4, iDn/10.);
                //cout << "Scales " << scaleUp << " "<< scaleDn << endl;

                GetDistancesPx(scaleUp, scaleDn);

                double dist = iterate(bestLayout);

                if(dist > minSepar)
                    goto gotoAfterScaleLoop;
        }
        
        gotoAfterScaleLoop:

        xx.resize(nLeg);
        yy.resize(nLeg);
        for(int i = 0; i < nLeg; ++i) {
            Point p;
            p.x = legBorders[i][bestLayout[i]].recs[0].fX;
            p.y = legBorders[i][bestLayout[i]].recs[0].fY;
            p = Px2NDC(p);
            xx[i] = p.x;
            yy[i] = p.y;
        }

        return {scaleUp, scaleDn};
        //cout << "done " << maxDist << endl;
    }




    void LoadHistoBorders(vector<Borders> &borders)
    {
        //Load histograms
        vector<TObject *> hists;
        for(const auto &&prim : *gPad->GetListOfPrimitives()) {
            cout << prim->GetName() <<" "<< prim->ClassName() <<  endl;
            if(strcmp(prim->GetName(),"hframe") && strcmp(prim->ClassName(),"TH1F"))
                if(dynamic_cast<TH1*>(prim)) 
                    hists.push_back(prim);
        }


        /////////////////////
        //GetDistance
        /////////////////////
        if(hists.size() == 0) {
            gPad->Modified();
            return;
        }

        for(unsigned i = 0; i < hists.size(); ++i) {
            Borders br;
            br.GetHistBorders(dynamic_cast<TH1*>(hists[i]));
            borders.push_back(br);
        }

    }


    double analyze(vector<int> &indx, vector<double> &dists) {
        double minDist = 1e40;
        //Distances to histograms
        for(unsigned k = 0; k < indx.size(); ++k) {
            minDist = min(minDist, distToHists[k][indx[k]]);
            dists[k] = distToHists[k][indx[k]];
            if(minDist < minSepar)
                return minDist;
        }

        for(unsigned k = 0;   k < indx.size(); ++k) 
        for(unsigned l = k+1; l < indx.size(); ++l) {
            double dist = MinDistance2(0.99*minSepar,
                           legBorders[k][indx[k]], legBorders[l][indx[l]]);
            minDist = min(minDist, dist);
            dists[k] = min(dists[k], dist);
            if(minDist < minSepar)
                return minDist;
        }


        return minDist;

    }


};

/// @name Auto-legend
/// Functions to place lagend automaticaly and prevent overlaps
///@{ 

/// Convert layout bit to string
///
/// Used to setup the position of particular TLegend by SetName
/// For example leg->SetName(SetLayout(kPos1)) for legend in to top left
/// Used internally by newLegend
/// @see newLegend()
/// @param pos the position bit, for instance kPos2 | kPos8
/// @return the number printed in decadic base to the string
inline const char *SetLayout(unsigned pos)
{
    return to_string(pos).c_str();
}


inline unsigned SimplifyPos(unsigned pos)
{
    if(pos & kPos5) return kPos5;
    if(pos & kPos2) pos &= !kPos1 & !kPos3;
    if(pos & kPos4) pos &= !kPos1 & !kPos7;
    if(pos & kPos6) pos &= !kPos3 & !kPos9;
    if(pos & kPos8) pos &= !kPos7 & !kPos9;
    return pos;
}

///
/// Struct containing NDC sizes and positions of the corresponding TLatex
///
struct RectangleNDC {
	double fX, fY, fWidth, fHeight;
	TLatex *lat;
	TText *tex;
	RectangleNDC() : lat(nullptr), tex(nullptr) {}
};


/// Calculates positions and height and width of rectangle which encapsulates TLatex.
///
/// Works well for 4 basic orientation, 0,90,180,270 degree.
/// Does not work for properly for vertical alignment to the bottom of the line
///
inline RectangleNDC GetNDC(TLatex *lat)
{

	int hAlign = lat->GetTextAlign() / 10 - 1;
	int vAlign = lat->GetTextAlign() % 10 - 1;
	if(vAlign == -1) vAlign = 0;

	RectangleNDC rTitle;

	rTitle.fHeight = lat->GetYsize() / (gPad->GetY2()-gPad->GetY1());
	rTitle.fWidth = lat->GetXsize()  / (gPad->GetX2()-gPad->GetX1());
	rTitle.fX = lat->GetX();
	rTitle.fY = lat->GetY();


	if(lat->GetTextAngle() == 0) {
		rTitle.fX -= hAlign*rTitle.fWidth/2.;
		rTitle.fY -= vAlign*rTitle.fHeight/2.;
	}
	else {
		swap(rTitle.fWidth, rTitle.fHeight);
		double rat = gPad->GetWh() / double(gPad->GetWw());
		rat /= gPad->GetAbsWNDC() / double(gPad->GetAbsHNDC());
		rTitle.fWidth *= rat;
		rTitle.fHeight *= 1./rat;

		if(lat->GetTextAngle() == 90) {
			rTitle.fY -= (hAlign)*rTitle.fHeight/2.;
			vAlign = 2-vAlign;
			rTitle.fX -=1.*vAlign*rTitle.fWidth/2.;
		}
		else if(lat->GetTextAngle() == 270) {
			hAlign = 2-hAlign;
			rTitle.fY -= (hAlign)*rTitle.fHeight/2.;
			rTitle.fX -=1.*vAlign*rTitle.fWidth/2.;
		}
	}
	rTitle.lat = lat;
	return rTitle;
}

/// Define new legend at position pos with nCols columns.
///
/// The text size is set to be the same as the size of y-axis title
/// All this values can be changed later on manually, like
/// leg->SetTextSize(1.3*leg->GetTextSize())
///
TLegend *newLegend(unsigned pos, int nCols = 1)
{
        TLegend *leg = new TLegend(0., 0., 0., 0.);
        leg->SetName(SetLayout(pos));
        leg->SetTextSize(GetYaxis()->GetTitleSize());
        leg->SetNColumns(nCols);
        return leg;
}

/// Calculates sizes of the legend in NDC
///
/// Note that due to the bug in root the y-size to be set to TLegend object
/// is sometimes diffrent than the real size
void GetLegendSizes(TLegend *leg, double &SizeX, double &SizeY, double &SizeYroot)
{

    double lineSeparation = 1.2;
    double markerSize = 1.5;
    double columnSeparation = 0;

    //double fontSize = RelFontToPx(GetYaxis()->GetTitleSize());
    double fontSize = RelFontToPx(leg->GetTextSize());

    gPad->Update();
    //cout <<"There is just start " <<  gPad->GetUymax() << " "<< gPad->GetUymin() << endl;

    leg->SetBorderSize(0);
    //leg->SetLineStyle(2);
    leg->SetFillStyle(0);
    leg->SetTextSize(PxFontToRel(fontSize));

    double fSizeNDCy = fontSize / (gPad->GetWh() * gPad->GetAbsHNDC());
    double fSizeNDCx = fontSize / (gPad->GetWw() * gPad->GetAbsWNDC());

    int nHeaders = 0;
    int nLegItems = 0;

    int nCols = leg->GetNColumns();
    vector<double> maxC(nCols, 0);
    double headerW = 0;
    int colId = 0;
    for(const auto &&entry : *leg->GetListOfPrimitives()) {
        TString lab = dynamic_cast<TLegendEntry *>(entry)->GetLabel();
        TString opt = dynamic_cast<TLegendEntry *>(entry)->GetOption();
        cout << "Label "<< leg->GetNRows()<<" " << lab << endl;
        TLatex *lat = new TLatex(0.5, 0.5, lab);
        lat->SetTextSize(leg->GetTextSize());
        lat->SetTextFont(leg->GetTextFont());

        double textW = GetNDC(lat).fWidth;
        if(opt == "h") {
            headerW = max(headerW, textW);
            ++nHeaders;
            colId = -1;
        }
        else {
            maxC[colId] = max(maxC[colId], textW);
            ++nLegItems;
        }
        colId = (colId != nCols-1) ? colId + 1 : 0;
        delete lat;
    }
    double maxW = accumulate(maxC.begin(),maxC.end(),0.);

    int nRows = nHeaders + (nLegItems+nCols-1)/nCols;
    //cout << "RADEK nRows " << nRows << endl;
    SizeY = lineSeparation*fSizeNDCy*nRows;

    headerW += 0.1*fSizeNDCx*markerSize;

    SizeX = max(nCols*fSizeNDCx*markerSize +
               (nCols-1)*fSizeNDCx*columnSeparation + maxW, headerW);
    leg->SetMargin( nCols*fSizeNDCx*markerSize / SizeX);

    leg->SetColumnSeparation((nCols-1)*fSizeNDCx*columnSeparation / SizeX);


    SizeYroot = lineSeparation*fSizeNDCy*leg->GetNRows();

}




/// Setup the position of legends provided in legs array
///
/// Note that the legends must be drawn manually
/// To draw legends directly, use DrawLegends
///
void PlaceLegends(vector<TLegend*> legs, bool keepRange = false)
{
    int nScaleSteps = keepRange ? 1 : 10;

    vector<double> sizesX, sizesY, sizesYroot;

    for(unsigned i = 0; i < legs.size(); ++i) {
        double SizeX, SizeY, SizeYroot;
        GetLegendSizes(legs[i], SizeX, SizeY, SizeYroot);
        sizesX.push_back(SizeX);
        sizesY.push_back(SizeY);
        sizesYroot.push_back(SizeYroot);
    }

    PLACER placer;

    placer.init(legs, sizesX, sizesY);

    vector<double> xx, yy;
    double scaleUp, scaleDn;
    tie(scaleUp, scaleDn) = placer.GetSolution(xx, yy, nScaleSteps);
    //return;
    
    for(unsigned i = 0; i < legs.size(); ++i) {
        legs[i]->SetX1(xx[i]);
        legs[i]->SetY1(yy[i] + sizesY[i] - sizesYroot[i]);
        legs[i]->SetX2(xx[i] + sizesX[i]);
        legs[i]->SetY2(yy[i] + sizesY[i]);
    }

    
    TH1 *frameNow = GetFrame();
    if(gPad->GetLogy()) {
        double Max = frameNow->GetMaximum();
        double Min = frameNow->GetMinimum();
        frameNow->SetMaximum(Min * pow(Max/Min, scaleUp));
        frameNow->SetMinimum(Max / pow(Max/Min, scaleDn));
    }
    else {
        double Max = frameNow->GetMaximum();
        double Min = frameNow->GetMinimum();
        frameNow->SetMaximum(Min + (Max-Min)*scaleUp);
        frameNow->SetMinimum(Max - (Max-Min)*scaleDn);
    }

    /*
    for(int i = 0; i < legs.size(); ++i) {
        TLine *l = new TLine;
        l->SetNDC();
        l->DrawLineNDC(xx[i], yy[i], xx[i]+sizesX[i], yy[i]);
        l->DrawLineNDC(xx[i], yy[i]+sizesY[i], xx[i]+sizesX[i], yy[i]+sizesY[i]);
        l->DrawLineNDC(xx[i], yy[i], xx[i], yy[i]+sizesY[i]);
        l->DrawLineNDC(xx[i]+sizesX[i], yy[i], xx[i]+sizesX[i], yy[i]+sizesY[i]);

    }
    */
}

/// Draw legends provided in vector legs in such a way that there is no overlap.
///
/// in case that keepRange=true the y-axis range is not touch even if the amount of space is not sufficiend.
/// The keepRange=true is very useful when a grid of histograms is plotted
///
void DrawLegends(vector<TLegend*> legs, bool keepRange=false)
{
    PlaceLegends(legs, keepRange);
    for(auto & leg : legs)
        leg->Draw();
}










static Point Abs2Px(Point p, double scaleUp, double scaleDn)
{
   	double valXmin = gPad->GetX1();
    double valXmax = gPad->GetX2();

    double valYmin = gPad->GetY2() - scaleDn*(gPad->GetY2()-gPad->GetY1());
    double valYmax = gPad->GetY1() + scaleUp*(gPad->GetY2()-gPad->GetY1());


	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();


    if(gPad->GetLogx()) p.x = log10(p.x);
    if(gPad->GetLogy()) p.y = log10(p.y);

    p.x =  (p.x - valXmin)/(valXmax-valXmin) * pxW;
    p.y =  (p.y - valYmin)/(valYmax-valYmin) * pxH;

    //cout << p.x <<" "<<  p.y << endl;
    return p;
}

static Point Px2NDC(Point p)
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();
    p.x = p.x/pxW;
    p.y = p.y/pxH;
    return p;
}



inline double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h)
{
    Borders bS;
    bS.recs.push_back( {x, y, w, h} );
    bS.FromNDC2px();

    double minDist = 1e40;
    for(Borders &b : bor) {
        double d = MinDistance2(0.99*minSkip, b, bS);
        minDist = min(minDist, d);
        if(minDist <= minSkip) return minDist;
    }
    return sqrt(minDist);
}


static double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip)
{
    double minDist = 1e40;
    for(Borders &b : bor) {
        double d = MinDistance2(0.99*minSkip, b, bSingle);
        minDist = min(minDist, d);
        if(minDist <= minSkip) return minDist;
    }
    return sqrt(minDist);
}



static double MinDistance2(double minSkip, const Borders &br1, const Borders &br2)
{
    double minDist = 1e40;
    for(const auto &r1 : br1.recs)
        for(const auto &r2 : br2.recs) {
            double d = Borders::Distance2(r1, r2);
            minDist = min(minDist, d);
            if(minDist <= minSkip) return minDist;
        }
    return minDist;

}

void Borders::FromNDC2px()
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();

    for(Rectangle &r : recs) { 
        r.fX *= pxW;
        r.fY *= pxH;
        r.fWidth  *= pxW;
        r.fHeight *= pxH;
    }
}

/*
void Borders::FromPx2NDC()
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();

    for(Rectangle &r : recs) { 
        r.fX /= pxW;
        r.fY /= pxH;
        r.fWidth  /= pxW;
        r.fHeight /= pxH;
    }
}
*/



void Borders::FromAbs2px(double scaleUp, double scaleDn)
{
    for(Rectangle &r : recs) { 
        Point p1(r.fX, r.fY);
        Point p2(r.fX+r.fWidth, r.fY+r.fHeight);

        p1 = Abs2Px(p1, scaleUp, scaleDn);
        p2 = Abs2Px(p2, scaleUp, scaleDn);
        r.fX = p1.x;
        r.fY = p1.y;
        r.fWidth  = p2.x-p1.x;
        r.fHeight = p2.y-p1.y;
    }
}

///
/// Return the borders of the provided histogram
///
void Borders::GetHistBorders(TH1 *h)
{
    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double xMin = h->GetBinLowEdge(i);
        double xMax = h->GetBinLowEdge(i) + h->GetBinWidth(i);
        double yMin = h->GetBinContent(i)-h->GetBinError(i);
        double yMax = h->GetBinContent(i)+h->GetBinError(i);

        const double lim  = 1e-30;

        if(gPad->GetLogx()) {
            xMin = max(lim, xMin);
            xMax = max(lim, xMax);
        }
        if(gPad->GetLogy()) {
            yMin = max(lim, yMin);
            yMax = max(lim, yMax);
        }
        recs.push_back({xMax, yMin, xMax-xMin, yMax-yMin});

        TString drawOpt = h->GetDrawOption();
        if(drawOpt.Contains("hist") && i > 1) {
            double yMinLeft = h->GetBinContent(i-1)-h->GetBinError(i-1);
            double yMaxLeft = h->GetBinContent(i-1)+h->GetBinError(i-1);
            if(gPad->GetLogy()) {
                yMinLeft = max(lim, yMinLeft);
                yMaxLeft = max(lim, yMaxLeft);
            }
            double yMinNow = min(yMin, yMinLeft);
            double yMaxNow = max(yMax, yMaxLeft);
            recs.push_back({xMin, yMinNow, 0., yMaxNow-yMinNow});
        }

    }
}

inline double hypot2(double x, double y) {return x*x+y*y;}

/// Distance between two rectangles, 0 if overlap
double Borders::Distance2(const Rectangle &r1, const Rectangle &r2)
{
    double x1  = r1.fX;
    double x1b = r1.fX + r1.fWidth;
    double x2  = r2.fX;
    double x2b = r2.fX + r2.fWidth;

    double y1  = r1.fY;
    double y1b = r1.fY + r1.fHeight;
    double y2  = r2.fY;
    double y2b = r2.fY + r2.fHeight;

    bool left = x2b < x1;
    bool right = x1b < x2;
    bool bottom = y2b < y1;
    bool top = y1b < y2;
    if(top && left)
        return hypot2(x2b-x1, y2-y1b);
    else if(left && bottom)
        return hypot2(x2b-x1, y2b-y1);
    else if(bottom && right)
        return hypot2(x2-x1b, y2b-y1);
    else if(right && top)
        return hypot2(x2-x1b, y2-y1b);
    else if(left)
        return hypot2(x1 - x2b,0.);
    else if(right)
        return hypot2(x2 - x1b,0.);
    else if(bottom)
        return hypot2(y1 - y2b,0.);
    else if(top)
        return hypot2(y2 - y1b,0.);
    else
        return 0.;
}

///@}


/// @name Miscellaneous
/// Utilities for redraw axes, divde pad and auto y-range calculation
///@{ 

/// Redraw the current frame and axis ticks.
///
/// Useful if the original ticks were covered during plotting
///
inline void UpdateFrame()
{
    gPad->Update();
    gPad->RedrawAxis();
    TFrame *fr = gPad->GetFrame();
    fr->SetFillStyle(0);
    fr->Draw();
    gPad->Update();
}


/// Calculate the automatic range of the vertical axis and apply it.
///
/// Useful when several histograms are plotted to the same frame. 
/// In such case the y-range is normally chosen according to the one which is plotted first (without "same")
/// This method select such a range which corresponds to union of all y-ranges  
/// When applying this method y-range is no more plotting order dependent
inline void CalcYaxisRange()
{
	TH1 *hFrame = GetFrame();
	TAxis *ax = hFrame->GetXaxis();

	int iFirst = ax->GetFirst();
	int iLast  = ax->GetLast();

	double xMin = ax->GetBinLowEdge(iFirst);
	double xMax = ax->GetBinLowEdge(iLast) + ax->GetBinWidth(iLast);


	//cout << "First " << hFrame->GetXaxis()->GetFirst() << endl;
	//cout << "Last  " << hFrame->GetXaxis()->GetLast() << endl;
	cout << xMin << " "<< xMax << endl;

	double widthMin,  widthMax;
	double centerMin, centerMax;
	//return ;
	double yMin = 1e40, yMax = -1e40;
	for(auto && obj : *gPad->GetListOfPrimitives()) {

		if(obj->InheritsFrom(TH1::Class()))
			((TH1*)obj)->GetXaxis()->SetRangeUser(xMin, xMax);
		else if(obj->InheritsFrom(TMultiGraph::Class()))
			((TMultiGraph*)obj)->GetXaxis()->SetRangeUser(xMin, xMax);
		else if(obj->InheritsFrom(TGraph::Class()))
			((TGraph*)obj)->GetXaxis()->SetRangeUser(xMin, xMax);
		else if(obj->InheritsFrom(THStack::Class()))
			((THStack*)obj)->GetXaxis()->SetRangeUser(xMin, xMax);


		if(obj->InheritsFrom(TH1::Class()) ||
			obj->InheritsFrom(TMultiGraph::Class()) ||
			obj->InheritsFrom(TGraph::Class())  ||
			obj->InheritsFrom(THStack::Class()) ) {

			obj->Paint();
			//double zoom = gPad->GetUymax() - gPad->GetUymin();
			double yMinNow = gPad->GetUymin();
			double yMaxNow = gPad->GetUymax();

			if(yMinNow < yMin) {
				yMin = yMinNow;
				widthMin = yMaxNow - yMinNow;
				centerMin = (yMaxNow + yMinNow)/2.;
			}
			if(yMaxNow > yMax) {
				yMax = yMaxNow;
				widthMax = yMaxNow - yMinNow;
				centerMax = (yMaxNow + yMinNow)/2.;
			}
		}
	}


	double width = yMax - yMin;


	double zoomMin = width / widthMin;
	double zoomMax = width / widthMax;
	
	if(yMin >= 0. && !gPad->GetLogy())
		yMin = max(0., centerMin - zoomMin * widthMin/2.);
	else
		yMin =  centerMin - zoomMin * widthMin/2.;
	yMax = centerMax + zoomMax * widthMax/2.;
	

	//cout <<"Minima "<< yMin <<" "<< yMax << endl;
	if(gPad->GetLogy()) {
		cout <<"PRAHA "<< yMin <<" "<< yMax << endl;
		hFrame->SetMinimum(pow(10,yMin));
		hFrame->SetMaximum(pow(10,yMax));
	}
	else {
		hFrame->SetMinimum(yMin);
		hFrame->SetMaximum(yMax);
	}
}


/// Construct the lattice of pads according to the provided x and y sizes
///
/// Comparing to classical Divide method this function left no empty space
/// between frames but left spaces for axes.
/// Note that without setting the offsets, font sizes and tics manually
/// the ticks differs between pads.
/// From this reason, we recommend calling SetFTO() in each pad
/// The margins of the original pad are kept
/// @param xDivs vector defining horizontal sizes of frames, the absolute normalisation does not matter
/// @param yDivs vector defining vertical sizes of frames, the absolute normalisation does not matter
inline void DividePad(vector<double> xDivs, vector<double> yDivs)
{
    double lMag = gPad->GetLeftMargin();
    double rMag = gPad->GetRightMargin();
    double tMag = gPad->GetTopMargin();
    double bMag = gPad->GetBottomMargin();

    int nx = xDivs.size();
    int ny = yDivs.size();

    double Nx = accumulate(xDivs.begin(), xDivs.end(), 0.);
    double Ny = accumulate(yDivs.begin(), yDivs.end(), 0.);
    vector<double> xPos(nx+1), yPos(ny+1);
    xPos[0] = yPos[0] = 0;
    for(int i = 1; i <= nx; ++i)
        xPos[i] = xPos[i-1] + xDivs[i-1]/Nx * (1-lMag-rMag);
    for(int i = 1; i <= ny; ++i)
        yPos[i] = yPos[i-1] + yDivs[i-1]/Ny * (1-tMag-bMag);

    
    for(int ix = 0; ix < nx; ++ix)
    for(int iy = 0; iy < ny; ++iy) {
        double xl = (ix == 0)    ? 0 : lMag+xPos[ix];
        double xr = (ix == nx-1) ? 1 : lMag+xPos[ix+1];

        double yt = (iy == 0)    ? 1 : 1-tMag-yPos[iy];
        double yb = (iy == ny-1) ? 0 : 1-tMag-yPos[iy+1];

        int id = iy*nx+ix+1;
        TPad *pad = new TPad(SF("%s_%d", gPad->GetName(), id), "", xl, yb, xr, yt);

        pad->SetLeftMargin( (ix == 0) ? lMag/(xPos[1]+lMag) : 0);
        pad->SetRightMargin( (ix == nx-1) ? rMag/( 1-lMag - xPos[nx-1]) : 0);
        pad->SetTopMargin( (iy == 0) ? tMag/(yPos[1]+tMag) : 0);
        pad->SetBottomMargin( (iy == ny-1) ? bMag/(1-tMag - yPos[ny-1]) : 0);

//        if(ix == nx-1) {
//            cout << "Hela : " <<  lMag + xPos[nx-1] << endl;
//            cout << "lMag rMag : " <<lMag <<" "<<  rMag << endl;
//            cout  << "Ratio " << rMag / (1-lMag-xPos[nx-1]) << endl;
//        }

        pad->SetNumber(id);
        pad->Draw();
    }
}


///@}



}

#endif
