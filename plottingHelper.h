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


Point Px2NDC(Point p);


const char *SetLayout(unsigned pos);

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
