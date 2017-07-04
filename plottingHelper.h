#ifndef __plottingHelper__
#define __plottingHelper__

#include "TCanvas.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TString.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TMarker.h"
#include "THStack.h"
#include "TMultiGraph.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <string>
#include <cassert>

/// The namespace of whole Plotting Helper utility
/// 
///
namespace PlottingHelper {

using namespace std;

#ifndef SF
    #define SF TString::Format 
#endif



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
    kPos9 = BIT(9), ///< Right bottom corner
    kPos2c= BIT(10),///< Center of top side
    kPos4c= BIT(11),///< Center of left side
    kPos6c= BIT(12),///< Center of right side
    kPos8c= BIT(13) ///< Center of bottom side
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

            if((pos & kPos2c)  && i==0 && j==last/2)
                goto after;
            if((pos & kPos4c)  && i==last/2 && j==0)
                goto after;
            if((pos & kPos6c)  && i==last/2 && j==last)
                goto after;
            if((pos & kPos8c)  && i==last && j==last/2)
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

TH1   *GetFrame();
TAxis *GetXaxis();
TAxis *GetYaxis();
TAxis *GetZaxis();

double PxFontToRel(double px);
double RelFontToPx(double rel);
double GetAxisFractionX();
double GetAxisFractionY();
double TickAbsToRelX(double tick);
double TickAbsToRelY(double tick);

void SetFonts(double pxX, double pxY = -1, double pxT = -1);
void SetTicks(double tickX, double tickY = -1);
void SetLabelOffsetX(double off);
void SetTitleOffsetX(double off);
void SetLabelOffsetY(double off);
void SetTitleOffsetY(double off);
void SetFontsTicks(double px, double tickX, double tickY = -1.);
void SetOffsets(double lX, double tX, double lY, double tY);
void SetFTO(vector<double> fonts, vector<double> ticks, vector<double> offsets);


void DrawLatex(TVirtualPad *pad1, TVirtualPad *pad2, double x, double y, TString text, double fSize=-1.0, TString style="");
void DrawLatex(TVirtualPad *pad, double x, double y, TString text, double fSize=-1.0, TString style="");
void DrawLatex(double x, double y, TString text, double fSize=-1.0, TString style="");
void DrawLatexLRTB(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize=-1.0, TString style="");

void DrawLatexUp(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexDown(TVirtualPad *pad1,TVirtualPad *pad2, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexRight(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexLeft(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize=-1.0, TString style="");

void DrawLatexUp(TVirtualPad *pad, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexDown(TVirtualPad *pad, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexRight(TVirtualPad *pad, double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexLeft(TVirtualPad *pad, double Offset, TString text, double fSize=-1.0, TString style="");

void DrawLatexUp(double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexDown(double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexRight(double Offset, TString text, double fSize=-1.0, TString style="");
void DrawLatexLeft(double Offset, TString text, double fSize=-1.0, TString style="");



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
double MinDistance2(double minSkip, const Borders &br1, const Borders &br2);

//Don't search for lower distance if minSkip reach
double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h);
double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip);


Point Px2NDC(Point p);

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

    void init(vector<TLegend *> legs, vector<double> &sizesX, vector<double> &sizesY);
    double iterate(vector<int> &bestLayout);
    void GetLegendsPositions();
    void GetDistancesPx(double scaleUp, double scaleDn);
    pair<double,double> GetSolution(vector<double> &xx, vector<double> &yy,
                                                         int nScaleSteps=10);
    void LoadHistoBorders(vector<Borders> &borders);
    double analyze(vector<int> &indx, vector<double> &dists);
};

const char *SetLayout(unsigned pos);
unsigned SimplifyPos(unsigned pos);

///
/// Struct containing NDC sizes and positions of the corresponding TLatex
///
struct RectangleNDC {
	double fX, fY, fWidth, fHeight;
	TLatex *lat;
	TText *tex;
	RectangleNDC() : lat(nullptr), tex(nullptr) {}
};


RectangleNDC GetNDC(TLatex *lat);
TLegend *newLegend(unsigned pos, int nCols = 1);
void GetLegendSizes(TLegend *leg, double &SizeX, double &SizeY, double &SizeYroot);
void PlaceLegends(vector<TLegend*> legs, bool keepRange = false);
void DrawLegends(vector<TLegend*> legs, bool keepRange=false);

Point Abs2Px(Point p, double scaleUp, double scaleDn);
Point Px2NDC(Point p);


double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h);
double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip);
double MinDistance2(double minSkip, const Borders &br1, const Borders &br2);

inline double hypot2(double x, double y) {return x*x+y*y;}


void UpdateFrame();
void CalcYaxisRange();
void SetLeftRight(double lMargin, double rMargin);
void SetTopBottom(double tMargin, double bMargin);
void DividePad(vector<double> xDivs, vector<double> yDivs);
void DivideTransparent(vector<double> divX, vector<double> divY, bool useMargins = true);
vector<double> merge(vector<double> v1, vector<double> v2={}, vector<double> v3={}, vector<double> v4={}, vector<double> v5={}, vector<double> v6={}, vector<double> v7={} );
vector<double> repeat(vector<double> x, int n);
vector<double> group(double frame, double space, int n);

}

#endif
