#include "plottingHelper.h"




/// The namespace of whole Plotting Helper utility
/// 
///
namespace PlottingHelper {

using namespace std;

/////////////////////////////////////////////
//Declaration of internal classes and structs
/////////////////////////////////////////////


// Struct specifying borders of some object which can be approximated by rectangles
// We use such approach to define position of histograms
struct Borders{
    vector<Rectangle> recs;
    void GetHistBorders(TH1 *h);
    void FromAbs2px(double scaleUp, double scaleDn);
    void FromNDC2px();
    static double Distance2(const Rectangle &r1, const Rectangle &r2);
};



//Contains internal methods for automatic placement of legends
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

// Iterator which goes over all positions specified in the constructor
//
// Note that bitwise "or" operator | can be used to select more complex layout
struct PosIterator {
    PosIterator(int _nSteps, unsigned Pos) {
        pos = Pos;
        nSteps = _nSteps;
        last = nSteps-1;
        iSave = -1; jSave = nSteps;
        //Iterate();
    }
    bool Iterate();

    int nSteps, last;
    int iSave, jSave;
    unsigned pos;
};




//Don't search for lower distance if minSkip reach
double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h);
double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip);

//Don't search for lower distance if minSkip reach
double MinDistance2(double minSkip, const Borders &br1, const Borders &br2);


inline double hypot2(double x, double y) {return x*x+y*y;}


/// @name Getters
/// Functions returning elements of the active frame
///@{ 

/// The method returns current frame on the active pad.
///
/// Useful for example to set maximum of the y-axis, GetFrame()->SetMaximum(5)
TH1 *GetFrame()
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
    if(!hobj) {cout << "No frame created in active pad" << endl;}

   return hobj;
}

/// Return x-axis of the active frame
TAxis *GetXaxis() {
    TH1 *frame = GetFrame();
    return frame ? frame->GetXaxis() : nullptr;
}

/// Return y-axis of the active frame
TAxis *GetYaxis() {
    TH1 *frame = GetFrame();
    return frame ? frame->GetYaxis() : nullptr;
}

/// Return z-axis of the active frame
TAxis *GetZaxis() {
    TH1 *frame = GetFrame();
    return frame ? frame->GetZaxis() : nullptr;
}  

///@} 



/// @name Units helpers
/// Functions to convert units and obtain axis fractions
///@{ 

/// Transform font size from px to the relative units
///
/// The returned size is dependent on the currently active pad - the gPad pointer
/// @param px The font size in px
/// @return The font size in relative units
double PxFontToRel(double px)
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
double RelFontToPx(double rel)
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();
	if(pxW < pxH) return rel*pxW;
	         else return rel*pxH;
}


/// Fraction of the pad width covered by x-axis
double GetAxisFractionX()
{
    //return (gPad->GetUxmax() - gPad->GetUxmin())/(gPad->GetX2() - gPad->GetX1());
    return (1. - gPad->GetLeftMargin() - gPad->GetRightMargin());
}

/// Fraction of the pad height covered by y-axis
double GetAxisFractionY()
{
    //return (gPad->GetUymax() - gPad->GetUymin())/(gPad->GetY2() - gPad->GetY1());
    return (1. - gPad->GetBottomMargin() - gPad->GetTopMargin());
}

/// Convert tick length in px of x-axis to the relative units
double TickAbsToRelX(double tick)
{
    double nom = GetAxisFractionX() * gPad->GetAbsHNDC() * gPad->GetWh();
	return tick/nom;
}


/// Convert tick length in px of y-axis to the relative units
double TickAbsToRelY(double tick)
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
void SetFonts(double pxX, double pxY, double pxT)
{
    TH1 *frame = GetFrame();
    if(!frame) return;

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
void SetTicks(double tickX, double tickY)
{
    TH1 *frame = GetFrame();
    if(!frame) return;

    //gPad->Update();
    if(tickY < 0) tickY = tickX;
	frame->GetXaxis()->SetTickSize(TickAbsToRelX(tickX) );
	frame->GetYaxis()->SetTickSize(TickAbsToRelY(tickY) );
}

/// Set offset of the x-axis labels.
///
/// The zero offset correspond to labels "standing" on top of the x-axis.
/// The offset is measured in units corresponds to height of numbers like "012..."
/// Therefore offset 1 correspond to labels touching x-axis by from the bottom.
void SetLabelOffsetX(double off)
{
    TH1 *frame = GetFrame();
    if(!frame) return;

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
void SetTitleOffsetX(double off)
{
    TH1 *frame = GetFrame();
    if(!frame) return;

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
void SetLabelOffsetY(double off)
{
    TH1 *frame = GetFrame();
    if(!frame) return;
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
void SetTitleOffsetY(double off)
{
    TH1 *frame = GetFrame();
    if(!frame) return;
    TAxis *ax = frame->GetYaxis();
    off *= 3/4.;//to vertical size of "0"
    ax->SetTitleOffset( off*  RelFontToPx(ax->GetTitleSize())/1.6 / ( gPad->GetWw()*ax->GetTitleSize()* gPad->GetAbsWNDC() ) );
}

/// Set of the fonts and ticks sizes 
///
/// Set size of x-,y-axis labels and titles to px (in pixels).
/// In addition the x and y tick lengths can be specified
void SetFontsTicks(double px, double tickX, double tickY)
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
void SetOffsets(double lX, double tX, double lY, double tY)
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
void SetFTO(vector<double> fonts, vector<double> ticks, vector<double> offsets)
{
    fonts.resize(3, -1);
    ticks.resize(2, -1);
    offsets.resize(4, 0);

    SetFonts(fonts[0], fonts[1], fonts[2]);
    SetTicks(ticks[0], ticks[1]);
    SetOffsets(offsets[0], offsets[1], offsets[2], offsets[3]);
}

///@}



/// @name Generalized titles
/// Functions to simplify drawing of latex captions related to the particular frame(s)
/// Especially useful in case of describing complex grid of frames.
/// For captions inside of frame consider also using of automatic legend.
/// Conventional axis titles can be easily "reproduced" with these methods
///@{ 

/// Draw latex in coordinates x,y give a hull of two frames in given TPads
///
/// The rectangular hull of two frames corresponding to pad1 and pad2 is used as an NDC-like coordinate system
/// For example the x-value starts from the most left side of frame1 and/or frame2 and ends in the most right side of these frames
/// Note that in current implementation given pads cannot contain any sub-pad = there must be only one frame in pad1 and one frame
/// in pad2.
/// @param pad1 The pad where the first frame is located
/// @param pad2 The pad where the second frame is located (can be the same as pad1)
/// @param x,y  Normalised coordinates between 0 and 1 defined by the hull of frame1 and frame2, 0,0 is on bottom left
/// @param text Latex which will be plotted, consider using #splitline #bf, #it, #font[], #color[].
/// @param fSize Font size in pixels, by default the font size is taken from the title of the first frame.
/// @param style The text alignment and orientation. The text-orientation is defined by rotating
/// of letter "v", that means "v,<,>,^". The default text alignment to the center horizontally and vertically.
/// The text can be also aligned to the left "l" or right "r", vertically to the bottom "b" or top "t"
void DrawLatex(TVirtualPad *pad1, TVirtualPad *pad2, double x, double y, TString text, double fSize, TString style)
{
    TVirtualPad *padOrg = gPad;
    assert(pad1->GetCanvas() == pad2->GetCanvas());
    TCanvas *can = pad1->GetCanvas();

    //Setting of style
    TLatex *tex = new TLatex();

    int hAlign = 2;
    int vAlign = 2;

    double angle;
    if(style.Contains('>')) {
        if(style.Contains('l')) vAlign = 1;
        if(style.Contains('r')) vAlign = 3;
        if(style.Contains('b')) hAlign = 1;
        if(style.Contains('t')) hAlign = 3;
        angle = 90;
    }
    else if(style.Contains('<')) {
        if(style.Contains('l')) vAlign = 3;
        if(style.Contains('r')) vAlign = 1;
        if(style.Contains('b')) hAlign = 3;
        if(style.Contains('t')) hAlign = 1;
        angle = 270;
    }
    else if(style.Contains('^')) {
        if(style.Contains('l')) hAlign = 3;
        if(style.Contains('r')) hAlign = 1;
        if(style.Contains('b')) vAlign = 3;
        if(style.Contains('t')) vAlign = 1;
        angle = 180;
    }
    else {
        if(style.Contains('l')) hAlign = 1;
        if(style.Contains('r')) hAlign = 3;
        if(style.Contains('b')) vAlign = 1;
        if(style.Contains('t')) vAlign = 3;
        angle = 0;
    }

    tex->SetTextAlign(10*hAlign + vAlign);
    tex->SetTextAngle(angle);
    tex->SetTextFont(42);

    double l[2], r[2], b[2], t[2];

    TVirtualPad *pads[] = {pad1, pad2};

    for(int i = 0; i < 2; ++i) {
        l[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*pads[i]->GetLeftMargin();
        r[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*(1-pads[i]->GetRightMargin());
        b[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*pads[i]->GetBottomMargin();
        t[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*(1-pads[i]->GetTopMargin());
    }

    double left  = min(l[0], l[1]);
    double right = max(r[0], r[1]);
    double bottom= min(b[0], b[1]);
    double top   = max(t[0], t[1]);

    double xGlob = left   + (right-left) *x;
    double yGlob = bottom + (top-bottom) *y;


    auto isBtw = [](double low, double x, double high) {return low < x && x < high;};
    bool isInside = false;
    for(int i = 0; (i < 2 && pads[0] != pads[1]) || (i<1); ++i) {
        if(isBtw(l[i],xGlob,r[i]) && isBtw(b[i],yGlob,t[i]) ) {
            pads[i]->cd();
            if(fSize < 0) {
                TH1 *frame = GetFrame();
                fSize = frame ? RelFontToPx(frame->GetTitleSize()) : 12;
            }
            tex->SetTextSize(PxFontToRel(fSize));

            double xFactor = (r[i] - l[i]) / (right - left);
            double yFactor = (t[i] - b[i]) / (top - bottom);
            assert(x / xFactor < 1);
            assert(y / yFactor < 1);


            double xLoc = pads[i]->GetLeftMargin() + x/xFactor *
                           (1 - pads[i]->GetLeftMargin() - pads[i]->GetRightMargin());
            double yLoc = pads[i]->GetBottomMargin() + y/yFactor *
                           (1 - pads[i]->GetBottomMargin() - pads[i]->GetTopMargin());

            //TLine *l = new TLine();
            //l->SetLineColor(kRed);
            //l->DrawLineNDC(xLoc, yLoc, xLoc + 0.1, yLoc);
            tex->DrawLatexNDC(xLoc, yLoc, text);
            isInside = true;
        }
    }

    if(!isInside) {
        if(fSize < 0) {
            pads[0]->cd();
            TH1 *frame = GetFrame();
            fSize = frame ? RelFontToPx(frame->GetTitleSize()) : 12;
        }
        can->cd();
        tex->SetTextSize(PxFontToRel(fSize));
        tex->DrawLatexNDC(xGlob, yGlob, text);
    }

    gPad->Update();
    padOrg->cd();
}


/// Draw latex at coordinates of the frame inside the corresponding TPad
///
/// For description see the DrawLatex method, here only the one pad is provided
void DrawLatex(TVirtualPad *pad, double x, double y, TString text, double fSize, TString style)
{
    DrawLatex(pad, pad, x, y, text, fSize, style);
}


/// Draw latex at coordinates of the frame corresponding to active TPad
///
/// The same as DrawLatex but used for active pad
void DrawLatex(double x, double y, TString text, double fSize, TString style)
{
    DrawLatex(gPad, x, y, text, fSize, style);
}


/// Draw latex with distant Offset to the border of the frame1 and frame2 hull
///
/// The offset coding is used to cover left, right, bottom and top scenarios
void DrawLatexLRTB(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize, TString style)
{
    TVirtualPad *padOrg = gPad;
    if(fSize < 0) {
        pad1->cd();
        TH1 *frame = GetFrame();
        fSize = frame ? RelFontToPx(frame->GetTitleSize()) : 12;
        padOrg->cd();
    }
    auto Close = [](double x, double val) {return abs(x-1000*val) <1000;};

    double x = 0.5, y = 0.5;


    double b[2], t[2], l[2], r[2];

    TVirtualPad *pads[] = {pad1, pad2};
    for(int i = 0; i < 2; ++i) {
        l[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*pads[i]->GetLeftMargin();
        r[i] = pads[i]->GetAbsXlowNDC() + pads[i]->GetAbsWNDC()*(1-pads[i]->GetRightMargin());
        b[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*pads[i]->GetBottomMargin();
        t[i] = pads[i]->GetAbsYlowNDC() + pads[i]->GetAbsHNDC()*(1-pads[i]->GetTopMargin());
    }

    double left  = min(l[0], l[1]);
    double right = max(r[0], r[1]);
    double bottom= min(b[0], b[1]);
    double top   = max(t[0], t[1]);

    if(Close(Offset,2) || Close(Offset, 8)) {
        if(style.Contains('l')) x = 0.0;
        if(style.Contains('r')) x = 1.0;
        double unit = (fSize * 3./4) / ( pad1->GetWh() * (top - bottom) );

        if(Close(Offset,2)) y = 1 + unit*(Offset - 2000);
        else                y = 0 - unit*(Offset - 8000);

    }
    else if(Close(Offset,4) || Close(Offset, 6)) {
        if(style.Contains('b')) y = 0.0;
        if(style.Contains('t')) y = 1.0;
        double unit = (fSize * 3./4) / ( pad1->GetWw() * (right - left) );

        if(Close(Offset,6)) x = 1 + unit*(Offset - 6000);
        else                x = 0 - unit*(Offset - 4000);
    }
    else
        assert(0);

    DrawLatex(pad1, pad2, x, y, text, fSize, style);
}

/// Draw latex up of frame1 and frame2
///
/// Draw the latex above the frames 1 and 2. The Offset is in units of text height of text. It is exactly the same unit
/// as used for the title offsets
void DrawLatexUp(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize, TString style) {
    DrawLatexLRTB(pad1, pad2, 2000 + Offset, text, fSize, style);
}

///Draw latex down of frame1 and frame2
void DrawLatexDown(TVirtualPad *pad1,TVirtualPad *pad2, double Offset, TString text, double fSize, TString style) {
    DrawLatexLRTB(pad1, pad2, 8000 + Offset, text, fSize, style);
}

///Draw latex right of frame1 and frame2
void DrawLatexRight(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize, TString style) {
    DrawLatexLRTB(pad1, pad2, 6000 + Offset, text, fSize, style);
}

///Draw latex left of frame1 and frame2
void DrawLatexLeft(TVirtualPad *pad1, TVirtualPad *pad2, double Offset, TString text, double fSize, TString style) {
    DrawLatexLRTB(pad1, pad2, 4000 + Offset, text, fSize, style);
}


///Draw latex up of frame corresponding to given pad
void DrawLatexUp(TVirtualPad *pad, double Offset, TString text, double fSize, TString style) {
    DrawLatexUp(pad, pad, Offset, text, fSize, style);
}

///Draw latex down of frame corresponding to given pad
void DrawLatexDown(TVirtualPad *pad, double Offset, TString text, double fSize, TString style) {
    DrawLatexDown(pad, pad, Offset, text, fSize, style);
}

///Draw latex right of frame corresponding to given pad
void DrawLatexRight(TVirtualPad *pad, double Offset, TString text, double fSize, TString style) {
    DrawLatexRight(pad, pad, Offset, text, fSize, style);
}

///Draw latex left of frame corresponding to given pad
void DrawLatexLeft(TVirtualPad *pad, double Offset, TString text, double fSize, TString style) {
    DrawLatexLeft(pad, pad, Offset, text, fSize, style);
}


///Draw latex up of frame corresponding to active pad
void DrawLatexUp(double Offset, TString text, double fSize, TString style) {
    DrawLatexUp(gPad, Offset, text, fSize, style);
}
///Draw latex down of frame corresponding to active pad
void DrawLatexDown(double Offset, TString text, double fSize, TString style) {
    DrawLatexDown(gPad, Offset, text, fSize, style);
}
///Draw latex right of frame corresponding to active pad
void DrawLatexRight(double Offset, TString text, double fSize, TString style) {
    DrawLatexRight(gPad, Offset, text, fSize, style);
}
///Draw latex left of frame corresponding to active pad
void DrawLatexLeft(double Offset, TString text, double fSize, TString style) {
    DrawLatexLeft(gPad, Offset, text, fSize, style);
}




///@}


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
const char *SetLayout(unsigned pos)
{
    return to_string(pos).c_str();
}


unsigned SimplifyPos(unsigned pos)
{
    if(pos & kPos5) return kPos5;
    if(pos & kPos2) pos &= !kPos1 & !kPos3;
    if(pos & kPos4) pos &= !kPos1 & !kPos7;
    if(pos & kPos6) pos &= !kPos3 & !kPos9;
    if(pos & kPos8) pos &= !kPos7 & !kPos9;
    return pos;
}


/// Calculates positions and height and width of rectangle which encapsulates TLatex.
///
/// Works well for 4 basic orientation, 0,90,180,270 degree.
/// Does not work for properly for vertical alignment to the bottom of the line
///
RectangleNDC GetNDC(TLatex *lat)
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
TLegend *newLegend(unsigned pos, int nCols)
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

    double lineSeparation = 1.0;
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
void PlaceLegends(vector<TLegend*> legs, bool keepRange)
{
    int nScaleSteps = keepRange ? 1 : 20;

    vector<double> sizesX, sizesY, sizesYroot;

    for(unsigned i = 0; i < legs.size(); ++i) {
        double SizeX, SizeY, SizeYroot;
        GetLegendSizes(legs[i], SizeX, SizeY, SizeYroot);
        sizesX.push_back(SizeX);
        sizesY.push_back(SizeY);
        sizesYroot.push_back(SizeYroot);
    }

    PLACER placer;
    gPad->Update();

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
        //double Max = frameNow->GetMaximum();
        //double Min = frameNow->GetMinimum();
        double Max = pow(10,gPad->GetUymax());
        double Min = pow(10,gPad->GetUymin());

        frameNow->SetMaximum(Min * pow(Max/Min, scaleUp));
        frameNow->SetMinimum(Max / pow(Max/Min, scaleDn));
    }
    else {
        double Max = gPad->GetUymax();
        double Min = gPad->GetUymin();

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
void DrawLegends(vector<TLegend*> legs, bool keepRange)
{
    PlaceLegends(legs, keepRange);
    for(auto & leg : legs)
        leg->Draw();
}


    bool PosIterator::Iterate() {
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






    //Initalize layOuts, dims, legBorders, borders
    void PLACER::init(vector<TLegend *> legs, vector<double> &sizesX, vector<double> &sizesY) {
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

    double PLACER::iterate(vector<int> &bestLayout) {
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


    void PLACER::GetLegendsPositions() {

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


    void PLACER::GetDistancesPx(double scaleUp, double scaleDn) {
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


    pair<double,double> PLACER::GetSolution(vector<double> &xx, vector<double> &yy,
                                                         int nScaleSteps) {
        
        double scaleUp=1., scaleDn=1.;

        vector<int> bestLayout;

        for(int sum = 0; sum < nScaleSteps; ++sum)
        for(int iUp = sum; iUp >= 0; --iUp) {
                int iDn = sum - iUp;
                if(!gPad->GetLogy() && iDn != 0) continue;
                scaleUp = pow(4, iUp/10.);
                scaleDn = pow(4, iDn/10.);

                GetDistancesPx(scaleUp, scaleDn);

                double dist = iterate(bestLayout);

                if(dist > minSepar)
                    goto gotoAfterScaleLoop;
        }
        cout << "No solution found for legend" << endl;
        
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




    void PLACER::LoadHistoBorders(vector<Borders> &borders)
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


    double PLACER::analyze(vector<int> &indx, vector<double> &dists) {
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



Point Abs2Px(Point p, double scaleUp, double scaleDn)
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

Point Px2NDC(Point p)
{
	double pxH = gPad->GetWh()*gPad->GetAbsHNDC();
	double pxW = gPad->GetWw()*gPad->GetAbsWNDC();
    p.x = p.x/pxW;
    p.y = p.y/pxH;
    return p;
}



double MinDistanceSingle(vector<Borders> &bor, double minSkip, double x, double y, double w, double h)
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


double MinDistanceSingle(vector<Borders> &bor, Borders bSingle, double minSkip)
{
    double minDist = 1e40;
    for(Borders &b : bor) {
        double d = MinDistance2(0.99*minSkip, b, bSingle);
        minDist = min(minDist, d);
        if(minDist <= minSkip) return minDist;
    }
    return sqrt(minDist);
}



double MinDistance2(double minSkip, const Borders &br1, const Borders &br2)
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
void UpdateFrame()
{
    gPad->Update();
    gPad->RedrawAxis();
    TFrame *fr = gPad->GetFrame();
    if(fr) {
        fr->SetFillStyle(0);
        fr->Draw();
    }
    else
        cout << "No frame at active pad" << endl;
    gPad->Update();
}


/// Calculate the automatic range of the vertical axis and apply it.
///
/// Useful when several histograms are plotted to the same frame. 
/// In such case the y-range is normally chosen according to the one which is plotted first (without "same")
/// This method select such a range which corresponds to union of all y-ranges  
/// When applying this method y-range is no more plotting order dependent
void CalcYaxisRange()
{
	TH1 *hFrame = GetFrame();
    if(!hFrame) return;
	TAxis *ax = hFrame->GetXaxis();

	int iFirst = ax->GetFirst();
	int iLast  = ax->GetLast();

	double xMin = ax->GetBinLowEdge(iFirst);
	double xMax = ax->GetBinLowEdge(iLast) + ax->GetBinWidth(iLast);


	double widthMin=0,  widthMax=0;
	double centerMin=0, centerMax=0;
	//return ;
	double yMin = 1e100, yMax = -1e100;
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
	

	if(gPad->GetLogy()) {
		hFrame->SetMinimum(pow(10,yMin));
		hFrame->SetMaximum(pow(10,yMax));
	}
	else {
		hFrame->SetMinimum(yMin);
		hFrame->SetMaximum(yMax);
	}
}

/// Set simultaneously left and right margin of active pad
///
/// Implemented to avoid the bug in ROOT which causes that
/// gPad->SetLeftMargin(0.5); gPad->SetRightMargin(0.1);
/// is not always the same as
/// gPad->SetRightMargin(0.1); gPad->SetLeftMargin(0.5);
void SetLeftRight(double lMargin, double rMargin)
{
    gPad->SetLeftMargin(0);
    gPad->SetRightMargin(0);
    gPad->SetLeftMargin(lMargin);
    gPad->SetRightMargin(rMargin);
}

/// Set simultaneously top and bottom margin of active pad
///
/// Implemented to avoid the bug in ROOT which causes that
/// gPad->SetTopMargin(0.5); gPad->SetBottomMargin(0.1);
/// is not always the same as
/// gPad->SetBottomMargin(0.1); gPad->SetTopMargin(0.5);
void SetTopBottom(double tMargin, double bMargin)
{
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0);
    gPad->SetTopMargin(tMargin);
    gPad->SetBottomMargin(bMargin);
}



////////////////////////////////////////////////////////////////////
/// Construct the lattice of pads according to the provided x and y sizes
///
/// Comparing to classical Divide method this function left no empty space
/// between frames but left spaces for axes.
/// Note that without setting the offsets, font sizes and tics manually
/// the ticks differs between pads.
/// From this reason, we recommend calling SetFTO() in each pad
/// The margins of the original pad are kept
/// This method can be nested is some of the created pads needs additional subdivision.
/// The resulting pads can be accessed by the standard cd() command as in ROOT Divide().
/// @param xDivs vector defining horizontal sizes of frames, the absolute normalisation does not matter
/// @param yDivs vector defining vertical sizes of frames, the absolute normalisation does not matter
///
void DividePad(vector<double> xDivs, vector<double> yDivs)
{
    TVirtualPad *orgPad = gPad;

    double lMag = orgPad->GetLeftMargin();
    double rMag = orgPad->GetRightMargin();
    double tMag = orgPad->GetTopMargin();
    double bMag = orgPad->GetBottomMargin();

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
        TPad *pad = new TPad(SF("%s_%d", orgPad->GetName(), id), "", xl, yb, xr, yt);

        double lNew =  (ix == 0) ? lMag/(xPos[1]+lMag) : 0;
        double rNew =  (ix == nx-1) ? rMag/( 1-lMag - xPos[nx-1]) : 0;

        double tNew =  (iy == 0) ? tMag/(yPos[1]+tMag) : 0;
        double bNew =   (iy == ny-1) ? bMag/(1-tMag - yPos[ny-1]) : 0;

        pad->cd();
        SetLeftRight(lNew, rNew);
        SetTopBottom(tNew, bNew);
        orgPad->cd();
        
//        if(ix == nx-1) {
//            cout << "Hela : " <<  lMag + xPos[nx-1] << endl;
//            cout << "lMag rMag : " <<lMag <<" "<<  rMag << endl;
//            cout  << "Ratio " << rMag / (1-lMag-xPos[nx-1]) << endl;
//        }

        pad->SetNumber(id);
        pad->Draw();
    }
}




////////////////////////////////////////////////////////////////////
/// Construct the lattice of frames according to the provided x and y sizes using transparency trick.
///
/// This function creates set of transparent overlapping pads, each of them covering whole canvas.
/// The margins of of of these pads are chosen in such a way that the lattice of frames is constructed.
/// Comparing to the DividePad method this method allows to include spaces between plotted frames.
/// In case that "zero" spaces are selected the frame structure would look the same only the axis
/// would not be cut off by the frame edge.
/// On the other hand one needs to ensure that nuisance exes are not plotted
/// (i.e. setting their font size to 0).
/// Note that argument vectors divX and divY must always contain odd number of elements (sizes).
/// Comparing to DividePad method this method can be run several time for the same pad.
/// In such a way more complex frame structure can be created than the regular lattice.
/// The resulting pads can be accessed by the standard cd() command as in ROOT Divide().
/// In case of more calls on the same pad the indexes continues where there previous ended.
///
/// @param xDivs vector defining horizontal sizes of frames and spaces between them,
///    the absolute normalisation does not matter. The structure is the following {space, frame, space, frame, space}
///    In case of useMargins=true the very left and right space is taken from pad margins.
/// @param yDivs the same as xDivs but in vertical coordinate
/// @param useMargins specify wherever the current pad margins should be used.
///
void DivideTransparent(vector<double> divX, vector<double> divY, bool useMargins)
{

    TVirtualPad *orgPad = gPad;

    if(divX.size() % 2 != 1 || int(divX.size()) < (3-2*useMargins) ) {
        cout << "Wrong divX= " << divX.size() << endl;
        return;
    }
    if(divY.size() % 2 != 1 || int(divY.size()) < (3-2*useMargins) ) {
        cout << "Wrong divY= " << divY.size() << endl;
        return;
    }



    double sumX = accumulate(divX.begin(), divX.end(), 0.0);
    double sumY = accumulate(divY.begin(), divY.end(), 0.0);

    if(useMargins == true) {
        double l = orgPad->GetLeftMargin();
        double r = orgPad->GetRightMargin();
        double t = orgPad->GetTopMargin();
        double b = orgPad->GetBottomMargin();
        

        vector<double> newX, newY;
        newX.push_back(l * sumX/(1-l-r));
        newX.insert(newX.end(), divX.begin(), divX.end());
        newX.push_back(r * sumX/(1-l-r));

        newY.push_back(t * sumY/(1-t-b));
        newY.insert(newY.end(), divY.begin(), divY.end());
        newY.push_back(b * sumY/(1-t-b));

        DivideTransparent(newX, newY, false);
        return;
    }



    vector<double> edgesX(divX.size()+1);
    vector<double> edgesY(divY.size()+1);

    edgesX[0] = edgesY[0] = 0;
    for(unsigned i = 1; i < edgesX.size(); ++i)
        edgesX[i] = edgesX[i-1] + divX[i-1]/sumX;
    for(unsigned i = 1; i < edgesY.size(); ++i)
        edgesY[i] = edgesY[i-1] + divY[i-1]/sumY;

    int nPadsX = (divX.size()-1)/2;
    //int nPadsY = (divY.size()-1)/2;

    int kStart = 1;
    while(orgPad->GetPad(kStart))
        ++kStart;
    --kStart;


    int i = 1;
    //for(int x = 1; x < edgesX.size()-1; x+=2)
    for(int y = 1; y < int(edgesY.size())-1; y+=2) 
    for(int x = edgesX.size()-3; x >= 1; x-=2) {
        //i = (nPadsY-1-(y-1)/2) * nPadsX + (nPadsX-1-(x-1)/2) + 1; 
        i = ((y-1)/2) * nPadsX + ((x-1)/2) + 1; 
        
        TPad *pad = new TPad(TString::Format("%s_%d", orgPad->GetName(), kStart+i), "", 0.0, 0.0, 1, 1);
        pad->SetFillStyle(0);

        
        double l = max(0.,edgesX[x]);
        double r = max(0.,1-edgesX[x+1]);
        double t = max(0.,edgesY[y]);
        double b = max(0.,1-edgesY[y+1]);

        //cout <<"LeftRight " << l <<" "<<  r <<" "<< (1-l-r) <<" "<< i<< endl;

        pad->cd();
        SetLeftRight(l, r);
        SetTopBottom(t, b);
        orgPad->cd();


        pad->SetNumber(kStart+i);
        pad->Draw();
        ++i;
    }
}

/// Merge vectors into one
///
/// In current version up to 7 vectors can be merged
vector<double> merge(vector<double> v1, vector<double> v2, vector<double> v3, vector<double> v4, vector<double> v5, vector<double> v6, vector<double> v7)
{
    vector<vector<double>> v = {v1, v2, v3, v4, v5, v6, v7};
    vector<double> res;
    for(unsigned i = 0; i < v.size(); ++i)
        res.insert(res.end(), v[i].begin(), v[i].end());
    return res;
}

/// Repeat the provided pattern
///
/// Returns the vector which includes n times the vector x
vector<double> repeat(vector<double> x, int n)
{
    vector<double> res = x;
    for(int i = 1; i < n; ++i)
        res.insert(res.end(), x.begin(), x.end());
    return res;
}

/// Construct group of frames
///
/// Returns vector describing layout of n identical frames which are separated by the same space size
/// The spaces before first space and after last frame are not included
vector<double> group(double frame, double space, int n)
{
    vector<double> res;
    for(int i = 0; i < n-1; ++i) {
        res.push_back(frame);
        res.push_back(space);
    }
    res.push_back(frame);
    return res;
}




///@}



}

//dummy function
void plottingHelper()
{
}
