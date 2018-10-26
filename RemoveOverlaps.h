#include "TLatex.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "THLimitsFinder.h"
#include "TTimeStamp.h"
#include "TF1.h"
#include <vector>
#include <algorithm>

static const Int_t kHori = BIT(9); //defined in TPad

struct RectangleNDC {
	double fX, fY, fWidth, fHeight;
	TLatex *lat;
	TText *tex;
	RectangleNDC() : lat(nullptr), tex(nullptr) {}
};


struct Point {
	Point(): x(0), y(0) {}
	Point(double _x, double _y): x(_x), y(_y) {}
	double x, y;
};


struct Edge {
	Point p1, p2;

  //Point dir;
	//double lenght;
};


inline void CopyLatexStyleNDC(TLatex *lat, TLatex *textaxis) {
	lat->SetNDC();
	lat->SetTextColor(textaxis->GetTextColor());
	lat->SetTextFont(textaxis->GetTextFont());
	lat->SetTextAlign(textaxis->GetTextAlign());
	lat->SetTextAngle(textaxis->GetTextAngle());
	lat->SetTextSize(textaxis->GetTextSize());
}
inline void CopyTextStyleNDC(TText *lat, TText *textaxis) {
	lat->SetNDC();
	lat->SetTextColor(textaxis->GetTextColor());
	lat->SetTextFont(textaxis->GetTextFont());
	lat->SetTextAlign(textaxis->GetTextAlign());
	lat->SetTextAngle(textaxis->GetTextAngle());
	lat->SetTextSize(textaxis->GetTextSize());
}


/*
RectangleNDC GetNDCmyold(TLatex *lat, const char * axis)
{
	int hAlign = lat->GetTextAlign() / 10 - 1;
	int vAlign = lat->GetTextAlign() % 10 - 1;
	if(vAlign == -1) vAlign = 0;


	RectangleNDC recNDC;
	recNDC.fHeight = lat->GetYsize() / (gPad->GetY2()-gPad->GetY1());
	recNDC.fWidth = lat->GetXsize()  / (gPad->GetX2()-gPad->GetX1());


	recNDC.fX = lat->GetX() - hAlign*recNDC.fWidth/2.;
	recNDC.fY = lat->GetY() - vAlign*recNDC.fHeight/2.;


	lat->SetName(TString::Format("LABEL%s", axis));
	recNDC.lat = lat;

	return recNDC;
}
*/

inline RectangleNDC GetNDCmy(TLatex *lat)
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
		std::swap(rTitle.fWidth, rTitle.fHeight);
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





class myTGaxis : public TGaxis {

public:
	std::vector<RectangleNDC> recsNDC;
	RectangleNDC rTitleNDC;
	RectangleNDC rExpNDC;
	RectangleNDC rTitleMyNDC;
	bool isVisual = true;

	double GetMinimalOffset();

void PaintAxis(double_t xmin, double_t ymin, double_t xmax, double_t ymax,
                       double_t &wmin, double_t &wmax, Int_t &ndiv,   Option_t *chopt="",
                       double_t gridlength=0, Bool_t drawGridOnly=kFALSE) override
{

   const char *where = "PaintAxis";

   double_t alfa, beta, ratio1, ratio2, grid_side;
   double_t axis_lengthN = 0;
   double_t axis_length0 = 0;
   double_t axis_length1 = 0;
   double_t axis_length;
   double_t atick[3];
   double_t tick_side;
   double_t charheight;
   double_t phil, phi, sinphi, cosphi, asinphi, acosphi;
   double_t binLow = 0.,  binLow2 = 0.,  binLow3 = 0.;
   double_t binHigh = 0., binHigh2 = 0., binHigh3 = 0.;
   double_t binWidth = 0., binWidth2 = 0., binWidth3 = 0.;
   double_t xpl1, xpl2, ypl1, ypl2;
   double_t xtick = 0;
   double_t xtick0, xtick1, dxtick=0;
   double_t ytick, ytick0, ytick1;
   double_t wlabel, dwlabel;
   double_t xfactor, yfactor;
   double_t xlabel, ylabel, dxlabel;
   double_t xone, xtwo;
   double_t rlab;
   double_t x0, x1, y0, y1, xx0, xx1, yy0, yy1;
   xx0 = xx1 = yy0 = yy1 = 0;
   double_t xxmin, xxmax, yymin, yymax;
   xxmin = xxmax = yymin = yymax = 0;
   double_t xlside, xmside;
   double_t ww, af, rne;
   double_t xx, yy;
   double_t xmnlog, x00, x11, h2, h2sav, axmul, y;
   Float_t chupxvsav, chupyvsav;
   double_t rtxw, rtyw;
   Int_t nlabels, nticks, nticks0, nticks1;
   Int_t i, j, k, l, decade, ltick;
   Int_t mside, lside;
   Int_t nexe  = 0;
   Int_t lnlen = 0;
   Int_t iexe, if1, if2, na, nf, ih1, ih2, nbinin, nch, kmod;
   Int_t optionLog,optionBlank,optionVert,optionPlus,optionMinus,optionUnlab,optionPara;
   Int_t optionDown,optionRight,optionLeft,optionCent,optionEqual,optionDecimals=0,optionDot;
   Int_t optionY,optionText,optionGrid,optionSize,optionNoopt,optionInt,optionM,optionUp,optionX;
   Int_t optionTime;
   Int_t first=0,last=0,labelnumber;
   Int_t xalign, yalign;
   Int_t nn1, nn2, nn3, n1a, n2a, n3a, nb2, nb3;
   Int_t nbins=10, n1aold, nn1old;
   n1aold = nn1old = 0;
   Int_t ndyn;
   Int_t nhilab = 0;
   Int_t idn;
   Bool_t flexe = 0;
   Bool_t flexpo,flexne;
   char *label;
   char *chtemp;
   char *coded;
   char chlabel[256];
   char kchtemp[256];
   char chcoded[8];
   TLine *linegrid;
   TString timeformat;
   TString typolabel;
   time_t timelabel;
   double_t timed, wTimeIni;
   struct tm* utctis;
   double_t rangeOffset = 0;

   double_t epsilon   = 1e-5;
   const double_t kPI = TMath::Pi();

   double_t rwmi = wmin;
   double_t rwma = wmax;
   chtemp = &kchtemp[0];
   label  = &chlabel[0];
   linegrid  = 0;

   fFunction = (TF1*)gROOT->GetFunction(fFunctionName.Data());

   Bool_t noExponent = TestBit(TAxis::kNoExponent);

// If moreLogLabels = kTRUE more Log Intermediate Labels are drawn.
   Bool_t moreLogLabels = TestBit(TAxis::kMoreLogLabels);

// the following parameters correspond to the pad range in NDC
// and the user's coordinates in the pad

   double_t padh   = gPad->GetWh()*gPad->GetAbsHNDC();
   double_t rwxmin = gPad->GetX1();
   double_t rwxmax = gPad->GetX2();
   double_t rwymin = gPad->GetY1();
   double_t rwymax = gPad->GetY2();

   if(strchr(chopt,'G')) optionLog  = 1;  else optionLog  = 0;
   if(strchr(chopt,'B')) optionBlank= 1;  else optionBlank= 0;
   if(strchr(chopt,'V')) optionVert = 1;  else optionVert = 0;
   if(strchr(chopt,'+')) optionPlus = 1;  else optionPlus = 0;
   if(strchr(chopt,'-')) optionMinus= 1;  else optionMinus= 0;
   if(strchr(chopt,'U')) optionUnlab= 1;  else optionUnlab= 0;
   if(strchr(chopt,'P')) optionPara = 1;  else optionPara = 0;
   if(strchr(chopt,'O')) optionDown = 1;  else optionDown = 0;
   if(strchr(chopt,'R')) optionRight= 1;  else optionRight= 0;
   if(strchr(chopt,'L')) optionLeft = 1;  else optionLeft = 0;
   if(strchr(chopt,'C')) optionCent = 1;  else optionCent = 0;
   if(strchr(chopt,'=')) optionEqual= 1;  else optionEqual= 0;
   if(strchr(chopt,'Y')) optionY    = 1;  else optionY    = 0;
   if(strchr(chopt,'T')) optionText = 1;  else optionText = 0;
   if(strchr(chopt,'W')) optionGrid = 1;  else optionGrid = 0;
   if(strchr(chopt,'S')) optionSize = 1;  else optionSize = 0;
   if(strchr(chopt,'N')) optionNoopt= 1;  else optionNoopt= 0;
   if(strchr(chopt,'I')) optionInt  = 1;  else optionInt  = 0;
   if(strchr(chopt,'M')) optionM    = 1;  else optionM    = 0;
   if(strchr(chopt,'0')) optionUp   = 1;  else optionUp   = 0;
   if(strchr(chopt,'X')) optionX    = 1;  else optionX    = 0;
   if(strchr(chopt,'t')) optionTime = 1;  else optionTime = 0;
   if(strchr(chopt,'.')) optionDot  = 1;  else optionDot  = 0;
   if (TestBit(TAxis::kTickPlus))     optionPlus  = 2;
   if (TestBit(TAxis::kTickMinus))    optionMinus = 2;
   if (TestBit(TAxis::kCenterLabels)) optionM     = 1;
   if (TestBit(TAxis::kDecimals))     optionDecimals = 1;
   if (!gStyle->GetStripDecimals())   optionDecimals = 1;
   if (fAxis) {
      if (fAxis->GetLabels()) {
         optionM    = 1;
         optionText = 1;
         ndiv = fAxis->GetLast()-fAxis->GetFirst()+1;
      }
      TList *ml = fAxis->GetModifiedLabels();
      if (ml) {
         fModLabs = ml;
         fNModLabs = fModLabs->GetSize();
      } else {
         fModLabs  = 0;
         fNModLabs = 0;
      }
   }
   if (ndiv < 0) {
      Error(where, "Invalid number of divisions: %d",ndiv);
      return;
   }

// Set the grid length

   if (optionGrid) {
      if (gridlength == 0) gridlength = 0.8;
      linegrid = new TLine();
      linegrid->SetLineColor(gStyle->GetGridColor());
      if (linegrid->GetLineColor() == 0) linegrid->SetLineColor(GetLineColor());
      linegrid->SetLineStyle(gStyle->GetGridStyle());
      linegrid->SetLineWidth(gStyle->GetGridWidth());
   }

// No labels if the axis label offset is big.
// In that case the labels are not visible anyway.

   if (GetLabelOffset() > 1.1 ) optionUnlab = 1;

// Determine time format

   Int_t idF = fTimeFormat.Index("%F");
   if (idF>=0) {
      timeformat = fTimeFormat(0,idF);
   } else {
      timeformat = fTimeFormat;
   }

   //GMT option
   if (fTimeFormat.Index("GMT")>=0) optionTime =2;

   // Determine the time offset and correct for time offset not being integer.
   double_t timeoffset =0;
   if (optionTime) {
      if (idF>=0) {
         Int_t lnF = fTimeFormat.Length();
         TString stringtimeoffset = fTimeFormat(idF+2,lnF);
         Int_t year, mm, dd, hh, mi, ss;
         if (sscanf(stringtimeoffset.Data(), "%d-%d-%d %d:%d:%d", &year, &mm, &dd, &hh, &mi, &ss) == 6) {
            //Get time offset in seconds since EPOCH:
            struct tm tp;
            tp.tm_year   = year-1900;
            tp.tm_mon    = mm-1;
            tp.tm_mday   = dd;
            tp.tm_hour   = hh;
            tp.tm_min    = mi;
            tp.tm_sec    = ss;
            tp.tm_isdst  = 0; //no DST for UTC (and forced to 0 in MktimeFromUTC function)
            timeoffset = TTimeStamp::MktimeFromUTC(&tp);

            // Add the time offset's decimal part if it is there
            Int_t ids   = stringtimeoffset.Index("s");
            if (ids >= 0) {
               Float_t dp;
               Int_t lns   = stringtimeoffset.Length();
               TString sdp = stringtimeoffset(ids+1,lns);
               sscanf(sdp.Data(),"%g",&dp);
               timeoffset += dp;
            }
         } else {
            Error(where, "Time offset has not the right format");
         }
      } else {
         timeoffset = gStyle->GetTimeOffset();
      }
      wmin += timeoffset - (int)(timeoffset);
      wmax += timeoffset - (int)(timeoffset);

      // correct for time offset at a good limit (min, hour, day, month, year)
      struct tm* tp0;
      time_t timetp = (time_t)((Long_t)(timeoffset));
      double_t range = wmax - wmin;
      Long_t rangeBase = 60;
      if (range>60)       rangeBase = 60*20;       // minutes
      if (range>3600)     rangeBase = 3600*20;     // hours
      if (range>86400)    rangeBase = 86400*20;    // days
      if (range>2419200)  rangeBase = 31556736;    // months (average # days)
      rangeOffset = (double_t) ((Long_t)(timeoffset)%rangeBase);
      if (range>31536000) {
         tp0 = gmtime(&timetp);
         tp0->tm_mon   = 0;
         tp0->tm_mday  = 1;
         tp0->tm_hour  = 0;
         tp0->tm_min   = 0;
         tp0->tm_sec   = 0;
         tp0->tm_isdst = 1; // daylight saving time is on.
         rangeBase = (timetp-mktime(tp0)); // years
         rangeOffset = (double_t) (rangeBase);
      }
      wmax += rangeOffset;
      wmin += rangeOffset;
   }

// Determine number of divisions 1, 2 and 3
   n1a   = ndiv%100;
   n2a   = (ndiv%10000 - n1a)/100;
   n3a   = ndiv/10000;
   nn3   = TMath::Max(n3a,1);
   nn2   = TMath::Max(n2a,1)*nn3;
   nn1   = TMath::Max(n1a,1)*nn2+1;
   nticks= nn1;

// Axis bining optimisation is ignored if:
//   - the first and the last label are equal
//   - the number of divisions is 0
//   - less than 1 primary division is requested
//   - logarithmic scale is requested

   if (wmin == wmax || ndiv == 0 || n1a <= 1 || optionLog) {
      optionNoopt = 1;
      optionInt   = 0;
   }

// Axis bining optimisation
   if ( (wmax-wmin) < 1 && optionInt) {
      Error(where, "option I not available");
      optionInt = 0;
   }
   if (!optionNoopt || optionInt ) {

// Primary divisions optimisation
// When integer labelling is required, Optimize is invoked first
// and only if the result is not an integer labelling, AdjustBinSize is invoked.

      THLimitsFinder::Optimize(wmin,wmax,n1a,binLow,binHigh,nbins,binWidth,fChopt.Data());
      if (optionInt) {
         if (binLow != double_t(int(binLow)) || binWidth != double_t(int(binWidth))) {
            AdjustBinSize(wmin,wmax,n1a,binLow,binHigh,nbins,binWidth);
         }
      }
      if ((wmin-binLow)  > epsilon) { binLow  += binWidth; nbins--; }
      if ((binHigh-wmax) > epsilon) { binHigh -= binWidth; nbins--; }
      if (xmax == xmin) {
         rtyw  = (ymax-ymin)/(wmax-wmin);
         xxmin = xmin;
         xxmax = xmax;
         yymin = rtyw*(binLow-wmin)  + ymin;
         yymax = rtyw*(binHigh-wmin) + ymin;
      }
      else {
         rtxw  = (xmax-xmin)/(wmax-wmin);
         xxmin = rtxw*(binLow-wmin)  + xmin;
         xxmax = rtxw*(binHigh-wmin) + xmin;
         if (ymax == ymin) {
            yymin = ymin;
            yymax = ymax;
         }
         else {
            alfa  = (ymax-ymin)/(xmax-xmin);
            beta  = (ymin*xmax-ymax*xmin)/(xmax-xmin);
            yymin = alfa*xxmin + beta;
            yymax = alfa*xxmax + beta;
         }
      }
      if (fFunction) {
         yymin = ymin;
         yymax = ymax;
         xxmin = xmin;
         xxmax = xmax;
      } else {
         wmin = binLow;
         wmax = binHigh;
      }

// Secondary divisions optimisation
      nb2 = n2a;
      if (!optionNoopt && n2a > 1 && binWidth > 0) {
         THLimitsFinder::Optimize(wmin,wmin+binWidth,n2a,binLow2,binHigh2,nb2,binWidth2,fChopt.Data());
      }

// Tertiary divisions optimisation
      nb3 = n3a;
      if (!optionNoopt && n3a > 1 && binWidth2 > 0) {
         THLimitsFinder::Optimize(binLow2,binLow2+binWidth2,n3a,binLow3,binHigh3,nb3,binWidth3,fChopt.Data());
      }
      n1aold = n1a;
      nn1old = nn1;
      n1a    = nbins;
      nn3    = TMath::Max(nb3,1);
      nn2    = TMath::Max(nb2,1)*nn3;
      nn1    = TMath::Max(n1a,1)*nn2+1;
      nticks = nn1;
   }

// Coordinates are normalized
	//std::cout << "My range " << rwxmin <<" " << rwxmax << std::endl;
	//std::cout << "My range " << xmin <<" " << xmax << std::endl;

   ratio1 = 1/(rwxmax-rwxmin);
   ratio2 = 1/(rwymax-rwymin);
   x0     = ratio1*(xmin-rwxmin);
   x1     = ratio1*(xmax-rwxmin);
   y0     = ratio2*(ymin-rwymin);
   y1     = ratio2*(ymax-rwymin);
   if (!optionNoopt || optionInt ) {
      xx0 = ratio1*(xxmin-rwxmin);
      xx1 = ratio1*(xxmax-rwxmin);
      yy0 = ratio2*(yymin-rwymin);
      yy1 = ratio2*(yymax-rwymin);
   }

   if ((x0 == x1) && (y0 == y1)) {
      Error(where, "length of axis is 0");
      return;
   }

// Return wmin, wmax and the number of primary divisions
   if (optionX) {
      ndiv = n1a;
      return;
   }

   Int_t maxDigits = fgMaxDigits;

   TLatex *textaxis = new TLatex();
   SetLineStyle(1); // axis line style
   textaxis->SetTextColor(GetTextColor());
   textaxis->SetTextFont(GetTextFont());

   if (!gPad->IsBatch()) {
      gVirtualX->GetCharacterUp(chupxvsav, chupyvsav);
      gVirtualX->SetClipOFF(gPad->GetCanvasID());
   }

// Compute length of axis
   axis_length = TMath::Sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
   if (axis_length == 0) {
      Error(where, "length of axis is 0");
      goto L210;
   }
   if (!optionNoopt || optionInt) {
      axis_lengthN = TMath::Sqrt((xx1-xx0)*(xx1-xx0)+(yy1-yy0)*(yy1-yy0));
      axis_length0 = TMath::Sqrt((xx0-x0)*(xx0-x0)+(yy0-y0)*(yy0-y0));
      axis_length1 = TMath::Sqrt((x1-xx1)*(x1-xx1)+(y1-yy1)*(y1-yy1));
      if (axis_lengthN < epsilon) {
         optionNoopt = 1;
         optionInt   = 0;
         wmin        = rwmi;
         wmax        = rwma;
         n1a         = n1aold;
         nn1         = nn1old;
         nticks      = nn1;
         if (optionTime) {
            wmin        += timeoffset - (int)(timeoffset) + rangeOffset;
            wmax        += timeoffset - (int)(timeoffset) + rangeOffset;
         }
      }
   }

   if (x0 == x1) {
      phi  = 0.5*kPI;
      phil = phi;
   } else {
            phi = TMath::ATan2((y1-y0),(x1-x0));
      Int_t px0 = gPad->UtoPixel(x0);
      Int_t py0 = gPad->VtoPixel(y0);
      Int_t px1 = gPad->UtoPixel(x1);
      Int_t py1 = gPad->VtoPixel(y1);
      if (x0 < x1) phil = TMath::ATan2(double_t(py0-py1), double_t(px1-px0));
      else         phil = TMath::ATan2(double_t(py1-py0), double_t(px0-px1));
   }
   cosphi  = TMath::Cos(phi);
   sinphi  = TMath::Sin(phi);
   acosphi = TMath::Abs(cosphi);
   asinphi = TMath::Abs(sinphi);
   if (acosphi <= epsilon) { acosphi = 0;  cosphi  = 0; }
   if (asinphi <= epsilon) { asinphi = 0;  sinphi  = 0; }

// mside positive, tick marks on positive side
// mside negative, tick marks on negative side
// mside zero, tick marks on both sides
// Default is positive except for vertical axis

   mside=1;
   if (x0 == x1 && y1 > y0)       mside = -1;
   if (optionPlus)                mside = 1;
   if (optionMinus)               mside = -1;
   if (optionPlus && optionMinus) mside = 0;
   xmside = mside;
   lside = -mside;
   if (optionEqual) lside = mside;
   if (optionPlus && optionMinus) {
      lside = -1;
      if (optionEqual) lside=1;
   }
   xlside = lside;

// Tick marks size
   if(xmside >= 0) tick_side = 1;
   else            tick_side = -1;
   if (optionSize) atick[0] = tick_side*axis_length*fTickSize;
   else            atick[0] = tick_side*axis_length*0.03;

   atick[1] = 0.5*atick[0];
   atick[2] = 0.5*atick[1];

// Set the side of the grid
   if ((x0 == x1) && (y1 > y0))  grid_side =-1;
   else                          grid_side = 1;

// Compute Values if Function is given
   if(fFunction) {
      rwmi = fFunction->Eval(wmin);
      rwma = fFunction->Eval(wmax);
      if(rwmi > rwma) {
         double_t t = rwma;
         rwma = rwmi;
         rwmi = t;
      }
   }

// Draw the axis if needed...
   if (!optionBlank) {
      xpl1 = x0;
      xpl2 = x1;
      ypl1 = y0;
      ypl2 = y1;
      PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
		//std::cout << "My life " << xpl1<<" "<< ypl1<<" "<< xpl2<<" "<< ypl2 <<  std::endl;
   }

// Draw axis title if it exists
   if (!drawGridOnly && strlen(GetTitle())) {
      textaxis->SetTextSize (GetTitleSize());
      charheight = GetTitleSize();
      if ((GetTextFont() % 10) > 2) {
         charheight = charheight/gPad->GetWh();
      }
      double_t toffset = GetTitleOffset();
//////if (toffset < 0.1) toffset = 1; // Negative offset should be allowed
      if (x1 == x0) ylabel = xlside*1.6*charheight*toffset;
      else          ylabel = xlside*1.3*charheight*toffset;
      if (y1 == y0) ylabel = xlside*1.6*charheight*toffset;
		double ylabelMy = 0.0*ylabel;//RADEK
      double_t axispos;
      if (TestBit(TAxis::kCenterTitle)) axispos = 0.5*axis_length;
      else                       axispos = axis_length;
      if (TestBit(TAxis::kRotateTitle)) {
         if (x1 >= x0) {
            if (TestBit(TAxis::kCenterTitle)) textaxis->SetTextAlign(22);
            else                              textaxis->SetTextAlign(12);
            Rotate(axispos,ylabel,cosphi,sinphi,x0,y0,xpl1,ypl1);
         } else {
            if (TestBit(TAxis::kCenterTitle)) textaxis->SetTextAlign(22);
         else                                 textaxis->SetTextAlign(32);
            Rotate(axispos,ylabel,cosphi,sinphi,x0,y0,xpl1,ypl1);
         }
         textaxis->PaintLatex(gPad->GetX1() + xpl1*(gPad->GetX2() - gPad->GetX1()),
                              gPad->GetY1() + ypl1*(gPad->GetY2() - gPad->GetY1()),
                              phil=(kPI+phil)*180/kPI,
                              GetTitleSize(),
                              GetTitle());




      } else {
			double xpl1My, ypl1My;//RADEK
         if (x1 >= x0) {
            if (TestBit(TAxis::kCenterTitle)) textaxis->SetTextAlign(22);
            else                              textaxis->SetTextAlign(32);
            Rotate(axispos,ylabel,cosphi,sinphi,x0,y0,xpl1,ypl1);
            Rotate(axispos,ylabelMy,cosphi,sinphi,x0,y0,xpl1My,ypl1My);//RADEK
         } else {
            if (TestBit(TAxis::kCenterTitle)) textaxis->SetTextAlign(22);
         else                                 textaxis->SetTextAlign(12);
            Rotate(axispos,ylabel,cosphi,sinphi,x0,y0,xpl1,ypl1);
            Rotate(axispos,ylabelMy,cosphi,sinphi,x0,y0,xpl1My,ypl1My);//RADEK
         }
			
         textaxis->PaintLatex(gPad->GetX1() + xpl1*(gPad->GetX2() - gPad->GetX1()),
                              gPad->GetY1() + ypl1*(gPad->GetY2() - gPad->GetY1()),
                              phil*180/kPI,
                              GetTitleSize(),
                              GetTitle());
			//Original title RADEK
			TLatex *lat = new TLatex(xpl1, ypl1, GetTitle());
			lat->SetNDC();
			lat->SetTextColor(textaxis->GetTextColor());
			lat->SetTextFont(textaxis->GetTextFont());
			lat->SetTextAlign(textaxis->GetTextAlign());
			lat->SetTextAngle(phil*180/kPI);
			lat->SetTextSize(GetTitleSize());
			lat->Draw();
			lat->SetName(TString("Title")+GetName());
			rTitleNDC = GetNDCmy(lat);

			//std::cout << "MONIKAAAAAAAAA " << lat->GetTextAlign() << std::endl;
			//std::cout << "KockaX " << rTitleNDC.fX << " "<< xpl1 << std::endl;
			//std::cout << "KockaY " << rTitleNDC.fY << " "<< ypl1 << std::endl;
			if(!isVisual) lat->Delete();

			//Shifted title RADEK
			TLatex *latMy = new TLatex(xpl1My, ypl1My, GetTitle());
			latMy->SetNDC();
			latMy->SetTextColor(textaxis->GetTextColor());
			latMy->SetTextFont(textaxis->GetTextFont());
			latMy->SetTextAlign(textaxis->GetTextAlign());
			latMy->SetTextAngle(phil*180/kPI);
			latMy->SetTextSize(GetTitleSize());
			latMy->Draw();
			latMy->SetName(TString("TitleMy")+GetName());
			rTitleMyNDC = GetNDCmy(latMy);
			//std::cout << "KockaX " << rTitleNDC.fX << " "<< xpl1 << std::endl;
			//std::cout << "KockaY " << rTitleNDC.fY << " "<< ypl1 << std::endl;
			latMy->Delete();
			//if(!isVisual) latMy->Delete();





      }
   }

// No bining

   if (ndiv == 0)goto L210;
   if (wmin == wmax) {
      Error(where, "wmin (%f) == wmax (%f)", wmin, wmax);
      goto L210;
   }

// Labels preparation:
// Get character height
// Compute the labels orientation in case of overlaps
// (with alphanumeric labels for horizontal axis).

   charheight = GetLabelSize();
   if (optionText && GetLabelFont()%10 != 3) charheight *= 0.66666;
   textaxis->SetTextFont(GetLabelFont());
   if ((GetLabelFont()%10 < 2) && optionLog) // force TLatex mode in PaintLatex
      textaxis->SetTextFont((Int_t)(GetLabelFont()/10)*10+2);
   textaxis->SetTextColor(GetLabelColor());
   textaxis->SetTextSize (charheight);
   textaxis->SetTextAngle(GetTextAngle());
   if (GetLabelFont()%10 > 2) {
      charheight /= padh;
   }
   if (!optionUp && !optionDown && !optionY && !optionUnlab) {
      if (!drawGridOnly && optionText && ((ymin == ymax) || (xmin == xmax))) {
         textaxis->SetTextAlign(32);
         optionText = 2;
         Int_t nl = fAxis->GetLast()-fAxis->GetFirst()+1;
         double_t angle     = 0;
         for (i=fAxis->GetFirst(); i<=fAxis->GetLast(); i++) {
            textaxis->SetText(0,0,fAxis->GetBinLabel(i));
            if (textaxis->GetXsize() < (xmax-xmin)/nl) continue;
            angle = -20;
            break;
         }
         for (i=fAxis->GetFirst(); i<=fAxis->GetLast(); i++) {
            if ((!strcmp(fAxis->GetName(),"xaxis") && !gPad->TestBit(kHori))
              ||(!strcmp(fAxis->GetName(),"yaxis") &&  gPad->TestBit(kHori))) {
               if (nl > 50) angle = 90;
               if (fAxis->TestBit(TAxis::kLabelsHori)) angle = 0;
               if (fAxis->TestBit(TAxis::kLabelsVert)) angle = 90;
               if (fAxis->TestBit(TAxis::kLabelsUp))   angle = 20;
               if (fAxis->TestBit(TAxis::kLabelsDown)) angle =-20;
               if (angle ==   0) textaxis->SetTextAlign(23);
               if (angle == -20) textaxis->SetTextAlign(12);
               double_t s = -3;
               if (ymin == gPad->GetUymax()) {
                  if (angle == 0) textaxis->SetTextAlign(21);
                  s = 3;
               }
					std::cout << "Here first" << std::endl;
               textaxis->PaintLatex(fAxis->GetBinCenter(i),
                                    ymin + s*fAxis->GetLabelOffset()*(gPad->GetUymax()-gPad->GetUymin()),
                                    angle,
                                    textaxis->GetTextSize(),
                                    fAxis->GetBinLabel(i));
            } else if ((!strcmp(fAxis->GetName(),"yaxis") && !gPad->TestBit(kHori))
                    || (!strcmp(fAxis->GetName(),"xaxis") &&  gPad->TestBit(kHori))) {
               double_t s = -3;
               if (xmin == gPad->GetUxmax()) {
                  textaxis->SetTextAlign(12);
                  s = 3;
               }
					std::cout << "Here second" << std::endl;
               textaxis->PaintLatex(xmin + s*fAxis->GetLabelOffset()*(gPad->GetUxmax()-gPad->GetUxmin()),
                                    fAxis->GetBinCenter(i),
                                    0,
                                    textaxis->GetTextSize(),
                                    fAxis->GetBinLabel(i));
            } else {
					std::cout << "Here third" << std::endl;
               textaxis->PaintLatex(xmin - 3*fAxis->GetLabelOffset()*(gPad->GetUxmax()-gPad->GetUxmin()),
                                    ymin +(i-0.5)*(ymax-ymin)/nl,
                                    0,
                                    textaxis->GetTextSize(),
                                    fAxis->GetBinLabel(i));
            }
         }
      }
   }

// Now determine orientation of labels on axis
   if (!gPad->IsBatch()) {
      if (cosphi > 0) gVirtualX->SetCharacterUp(-sinphi,cosphi);
      else            gVirtualX->SetCharacterUp(sinphi,-cosphi);
      if (x0 == x1)   gVirtualX->SetCharacterUp(0,1);
      if (optionVert) gVirtualX->SetCharacterUp(0,1);
      if (optionPara) gVirtualX->SetCharacterUp(-sinphi,cosphi);
      if (optionDown) gVirtualX->SetCharacterUp(cosphi,sinphi);
   }

// Now determine text alignment
   xalign = 2;
   yalign = 1;
   if (x0 == x1)    xalign = 3;
   if (y0 != y1)    yalign = 2;
   if (optionCent)  xalign = 2;
   if (optionRight) xalign = 3;
   if (optionLeft)  xalign = 1;
   if (TMath::Abs(cosphi) > 0.9) {
      xalign = 2;
   } else {
      if (cosphi*sinphi > 0)  xalign = 1;
      if (cosphi*sinphi < 0)  xalign = 3;
   }
   textaxis->SetTextAlign(10*xalign+yalign);

// Position of labels in Y
   if (x0 == x1) {
      if (optionPlus && !optionMinus) {
         if (optionEqual) ylabel =  fLabelOffset/2 + atick[0];
         else             ylabel = -fLabelOffset;
      } else {
         ylabel = fLabelOffset;
         if (lside < 0)  ylabel += atick[0];
      }
   } else if (y0 == y1) {
      if (optionMinus && !optionPlus) {
         if ((GetLabelFont() % 10) == 3 ) {
            ylabel = fLabelOffset+0.5*
            ((gPad->AbsPixeltoY(0)-gPad->AbsPixeltoY((Int_t)fLabelSize))/
            (gPad->GetY2() - gPad->GetY1()));
         } else {
            ylabel = fLabelOffset+0.5*fLabelSize;
         }
         ylabel += TMath::Abs(atick[0]);
      } else {
         ylabel = -fLabelOffset;
         if (mside <= 0) ylabel -= TMath::Abs(atick[0]);
      }
      if (optionLog)  ylabel -= 0.5*charheight;
   } else {
      if (mside+lside >= 0) ylabel =  fLabelOffset;
      else                  ylabel = -fLabelOffset;
   }
   if (optionText) ylabel /= 2;

// Draw the linear tick marks if needed...
   if (!optionLog) {
      if (ndiv) {
         if (fFunction) {
            dxtick=(binHigh-binLow)/double_t(nticks-1);
         } else {
            if (optionNoopt && !optionInt) dxtick=axis_length/double_t(nticks-1);
            else                           dxtick=axis_lengthN/double_t(nticks-1);
         }
         for (k=0;k<nticks; k++) {
            ltick = 2;
            if (k%nn3 == 0) ltick = 1;
            if (k%nn2 == 0) ltick = 0;
            if (fFunction) {
               double_t xf = binLow+double_t(k)*dxtick;
               double_t zz = fFunction->Eval(xf)-rwmi;
               xtick = zz* axis_length / TMath::Abs(rwma-rwmi);
            } else {
               xtick = double_t(k)*dxtick;
            }
            ytick = 0;
            if (!mside) ytick -= atick[ltick];
            if ( optionNoopt && !optionInt) {
               Rotate(xtick,ytick,cosphi,sinphi,x0,y0,xpl2,ypl2);
               Rotate(xtick,atick[ltick],cosphi,sinphi,x0,y0,xpl1,ypl1);
            }
            else {
               Rotate(xtick,ytick,cosphi,sinphi,xx0,yy0,xpl2,ypl2);
               Rotate(xtick,atick[ltick],cosphi,sinphi,xx0,yy0,xpl1,ypl1);
            }
            if (optionVert) {
               if ((x0 != x1) && (y0 != y1)) {
                  if (mside) {
                     xpl1 = xpl2;
                     if (cosphi > 0) ypl1 = ypl2 + atick[ltick];
                     else            ypl1 = ypl2 - atick[ltick];
                  }
                  else {
                     xpl1 = 0.5*(xpl1 + xpl2);
                     xpl2 = xpl1;
                     ypl1 = 0.5*(ypl1 + ypl2) + atick[ltick];
                     ypl2 = 0.5*(ypl1 + ypl2) - atick[ltick];
                  }
               }
            }
            if (!drawGridOnly) PaintLineNDC(xpl1, ypl1, xpl2, ypl2);

            if (optionGrid) {
               if (ltick == 0) {
                  if (optionNoopt && !optionInt) {
                     Rotate(xtick,0,cosphi,sinphi,x0,y0 ,xpl2,ypl2);
                     Rotate(xtick,grid_side*gridlength ,cosphi,sinphi,x0,y0 ,xpl1,ypl1);
                  }
                  else {
                     Rotate(xtick,0,cosphi ,sinphi,xx0,yy0 ,xpl2,ypl2);
                     Rotate(xtick,grid_side*gridlength ,cosphi,sinphi,xx0,yy0 ,xpl1,ypl1);
                  }
                  linegrid->PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
               }
            }
         }
         xtick0 = 0;
         xtick1 = xtick;

         if (fFunction) axis_length0 = binLow-wmin;
         if ((!optionNoopt || optionInt) && axis_length0) {
            nticks0 = Int_t(axis_length0/dxtick);
            if (nticks0 > 1000) nticks0 = 1000;
            for (k=0; k<=nticks0; k++) {
               ltick = 2;
               if (k%nn3 == 0) ltick = 1;
               if (k%nn2 == 0) ltick = 0;
               ytick0 = 0;
               if (!mside) ytick0 -= atick[ltick];
               if (fFunction) {
                  xtick0 = (fFunction->Eval(binLow - double_t(k)*dxtick)-rwmi)
                           * axis_length / TMath::Abs(rwma-rwmi);
               }
               Rotate(xtick0,ytick0,cosphi,sinphi,xx0,yy0 ,xpl2,ypl2);
               Rotate(xtick0,atick[ltick],cosphi,sinphi,xx0,yy0 ,xpl1,ypl1);
               if (optionVert) {
                  if ((x0 != x1) && (y0 != y1)) {
                     if (mside) {
                        xpl1 = xpl2;
                        if (cosphi > 0) ypl1 = ypl2 + atick[ltick];
                        else            ypl1 = ypl2 - atick[ltick];
                     }
                     else {
                        xpl1 = 0.5*(xpl1 + xpl2);
                        xpl2 = xpl1;
                        ypl1 = 0.5*(ypl1 + ypl2) + atick[ltick];
                        ypl2 = 0.5*(ypl1 + ypl2) - atick[ltick];
                     }
                  }
               }
               if (!drawGridOnly) PaintLineNDC(xpl1, ypl1, xpl2, ypl2);

               if (optionGrid) {
                  if (ltick == 0) {
                     Rotate(xtick0,0,cosphi,sinphi,xx0,yy0,xpl2,ypl2);
                     Rotate(xtick0,grid_side*gridlength ,cosphi,sinphi,xx0,yy0 ,xpl1,ypl1);
                     linegrid->PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
                  }
               }
               xtick0 -= dxtick;
            }
         }

         if (fFunction) axis_length1 = wmax-binHigh;
         if ((!optionNoopt || optionInt) && axis_length1) {
            nticks1 = int(axis_length1/dxtick);
            if (nticks1 > 1000) nticks1 = 1000;
            for (k=0; k<=nticks1; k++) {
               ltick = 2;
               if (k%nn3 == 0) ltick = 1;
               if (k%nn2 == 0) ltick = 0;
               ytick1 = 0;
               if (!mside) ytick1 -= atick[ltick];
               if (fFunction) {
                  xtick1 = (fFunction->Eval(binHigh + double_t(k)*dxtick)-rwmi)
                           * axis_length / TMath::Abs(rwma-rwmi);
               }
               Rotate(xtick1,ytick1,cosphi,sinphi,xx0,yy0 ,xpl2,ypl2);
               Rotate(xtick1,atick[ltick],cosphi,sinphi,xx0,yy0 ,xpl1,ypl1);
               if (optionVert) {
                  if ((x0 != x1) && (y0 != y1)) {
                     if (mside) {
                        xpl1 = xpl2;
                        if (cosphi > 0) ypl1 = ypl2 + atick[ltick];
                        else            ypl1 = ypl2 - atick[ltick];
                     }
                     else {
                        xpl1 = 0.5*(xpl1 + xpl2);
                        xpl2 = xpl1;
                        ypl1 = 0.5*(ypl1 + ypl2) + atick[ltick];
                        ypl2 = 0.5*(ypl1 + ypl2) - atick[ltick];
                     }
                  }
               }
               if (!drawGridOnly) PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
               if (optionGrid) {
                  if (ltick == 0) {
                     Rotate(xtick1,0,cosphi,sinphi,xx0,yy0 ,xpl2,ypl2);
                     Rotate(xtick1,grid_side*gridlength,cosphi,sinphi,xx0,yy0,xpl1,ypl1);
                     linegrid->PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
                  }
               }
               xtick1 += dxtick;
            }
         }
      }
   }

// Draw the numeric labels if needed...
   if (!drawGridOnly && !optionUnlab) {
      if (!optionLog) {
         if (n1a) {
// Spacing of labels
            if ((wmin == wmax) || (ndiv == 0)) {
               Error(where, "wmin (%f) == wmax (%f), or ndiv == 0", wmin, wmax);
               goto L210;
            }
            wlabel  = wmin;
            dwlabel = (wmax-wmin)/double_t(n1a);
            if (optionNoopt && !optionInt) dxlabel = axis_length/double_t(n1a);
            else                           dxlabel = axis_lengthN/double_t(n1a);

            if (!optionText && !optionTime) {

// We have to decide what format to generate
// (for numeric labels only)
// Test the magnitude, decide format
               flexe  = kFALSE;
               nexe   = 0;
               flexpo = kFALSE;
               flexne = kFALSE;
               ww     = TMath::Max(TMath::Abs(wmin),TMath::Abs(wmax));

// First case : (wmax-wmin)/n1a less than 0.001
// (0.001 fgMaxDigits of 5 (fgMaxDigits) characters). Then we use x 10 n
// format. If af >=0 x10 n cannot be used
               double_t xmicros = 0.00099;
               if (maxDigits) xmicros = TMath::Power(10,-maxDigits);
               if (!noExponent && (TMath::Abs(wmax-wmin)/double_t(n1a)) < xmicros) {
                  af    = TMath::Log10(ww) + epsilon;
                  if (af < 0) {
                     flexe   = kTRUE;
                     nexe    = int(af);
                     iexe    = TMath::Abs(nexe);
                     if (iexe%3 == 1)     iexe += 2;
                     else if(iexe%3 == 2) iexe += 1;
                     if (nexe < 0) nexe = -iexe;
                     else          nexe =  iexe;
                     wlabel  = wlabel*TMath::Power(10,iexe);
                     dwlabel = dwlabel*TMath::Power(10,iexe);
                     if1     = maxDigits;
                     if2     = maxDigits-2;
                     goto L110;
                  }
               }
               if (ww >= 1) af = TMath::Log10(ww);
               else         af = TMath::Log10(ww*0.0001);
               af += epsilon;
               nf  = Int_t(af)+1;
               if (!noExponent && nf > maxDigits)  flexpo = kTRUE;
               if (!noExponent && nf < -maxDigits) flexne = kTRUE;

// Use x 10 n format. (only powers of 3 allowed)

               if (flexpo) {
                  flexe = kTRUE;
                  while (1) {
                     nexe++;
                     ww      /= 10;
                     wlabel  /= 10;
                     dwlabel /= 10;
                     if (nexe%3 == 0 && ww <= TMath::Power(10,maxDigits-1)) break;
                  }
               }

               if (flexne) {
                  flexe = kTRUE;
                  rne   = 1/TMath::Power(10,maxDigits-2);
                  while (1) {
                     nexe--;
                     ww      *= 10;
                     wlabel  *= 10;
                     dwlabel *= 10;
                     if (nexe%3 == 0 && ww >= rne) break;
                  }
               }

               na = 0;
               for (i=maxDigits-1; i>0; i--) {
                  if (TMath::Abs(ww) < TMath::Power(10,i)) na = maxDigits-i;
               }
               ndyn = n1a;
               while (ndyn) {
                  double_t wdyn = TMath::Abs((wmax-wmin)/ndyn);
                  if (wdyn <= 0.999 && na < maxDigits-2) {
                     na++;
                     ndyn /= 10;
                  }
                  else break;
               }

               if2 = na;
               if1 = TMath::Max(nf+na,maxDigits)+1;
L110:
               if (TMath::Min(wmin,wmax) < 0)if1 = if1+1;
               if1 = TMath::Min(if1,32);

// In some cases, if1 and if2 are too small....
               while (dwlabel < TMath::Power(10,-if2)) {
                  if1++;
                  if2++;
               }
               coded = &chcoded[0];
               if (if1 > 14) if1=14;
               if (if2 > 14) if2=14;
               if (if2>0) snprintf(coded,8,"%%%d.%df",if1,if2);
               else       snprintf(coded,8,"%%%d.%df",if1+1,1);
            }

// We draw labels

            snprintf(chtemp,256,"%g",dwlabel);
            Int_t ndecimals = 0;
            if (optionDecimals) {
               char *dot = strchr(chtemp,'.');
               if (dot) {
                  ndecimals = chtemp + strlen(chtemp) -dot;
               } else {
                  char *exp;
                  exp = strstr(chtemp,"e-");
                  if (exp) {
                     sscanf(&exp[2],"%d",&ndecimals);
                     ndecimals++;
                  }
               }
            }
            if (optionM) nlabels = n1a-1;
            else         nlabels = n1a;
            wTimeIni = wlabel;
            for ( k=0; k<=nlabels; k++) {
               if (fFunction) {
                  double_t xf = binLow+double_t(k*nn2)*dxtick;
                  double_t zz = fFunction->Eval(xf)-rwmi;
                  wlabel = xf;
                  xlabel = zz* axis_length / TMath::Abs(rwma-rwmi);
               } else {
                  xlabel = dxlabel*k;
               }
               if (optionM)    xlabel += 0.5*dxlabel;

               if (!optionText && !optionTime) {
                  snprintf(label,256,&chcoded[0],wlabel);
                  label[28] = 0;
                  wlabel += dwlabel;

                  LabelsLimits(label,first,last);  //Eliminate blanks

                  if (label[first] == '.') { //check if '.' is preceded by a digit
                     strncpy(chtemp, "0",256);
                     strlcat(chtemp, &label[first],256);
                     strncpy(label, chtemp,256);
                     first = 1; last = strlen(label);
                  }
                  if (label[first] == '-' && label[first+1] == '.') {
                     strncpy(chtemp, "-0",256);
                     strlcat(chtemp, &label[first+1],256);
                     strncpy(label, chtemp, 256);
                     first = 1; last = strlen(label);
                  }

// We eliminate the non significant 0 after '.'
                  if (ndecimals) {
                     char *adot = strchr(label,'.');
                     if (adot) adot[ndecimals] = 0;
                  } else {
                     while (label[last] == '0') { label[last] = 0; last--;}
                  }

// We eliminate the dot, unless dot is forced.
                  if (label[last] == '.') {
                     if (!optionDot) { label[last] = 0; last--;}
                  }

// Make sure the label is not "-0"
                  if (last-first == 1 && label[first] == '-'
                                      && label[last]  == '0') {
                     strncpy(label, "0", 256);
                     label[last] = 0;
                  }
               }

// Generate the time labels

               if (optionTime) {
                  timed = wlabel + (int)(timeoffset) - rangeOffset;
                  timelabel = (time_t)((Long_t)(timed));
                  if (optionTime == 1) {
                     utctis = localtime(&timelabel);
                  } else {
                     utctis = gmtime(&timelabel);
                  }
                  TString timeformattmp;
                  if (timeformat.Length() < 220) timeformattmp = timeformat;
                  else timeformattmp = "#splitline{Format}{too long}";

// Appends fractional part if seconds displayed
                  if (dwlabel<0.9) {
                     double tmpdb;
                     int tmplast;
                     snprintf(label, 256, "%%S%7.5f", modf(timed,&tmpdb));
                     tmplast = strlen(label)-1;

// We eliminate the non significant 0 after '.'
                     while (label[tmplast] == '0') {
                        label[tmplast] = 0; tmplast--;
                     }

                     timeformattmp.ReplaceAll("%S",label);
// replace the "0." at the beginning by "s"
                     timeformattmp.ReplaceAll("%S0.","%Ss");

                  }

                  if (utctis != nullptr) {
                     strftime(label, 256, timeformattmp.Data(), utctis);
                  } else {
                     strncpy(label, "invalid", 256);
                  }
                  strncpy(chtemp, &label[0], 256);
                  first = 0; last=strlen(label)-1;
                  wlabel = wTimeIni + (k+1)*dwlabel;
               }

// We generate labels (numeric or alphanumeric).

               if (optionNoopt && !optionInt)
                        Rotate (xlabel,ylabel,cosphi,sinphi,x0,y0,xx,yy);
               else     Rotate (xlabel,ylabel,cosphi,sinphi,xx0,yy0,xx,yy);
               if (y0 == y1 && !optionDown && !optionUp) {
                  yy -= 0.80*charheight;
               }
               if (optionVert) {
                  if (x0 != x1 && y0 != y1) {
                     if (optionNoopt && !optionInt)
                           Rotate (xlabel,0,cosphi,sinphi,x0,y0,xx,yy);
                     else  Rotate (xlabel,0,cosphi,sinphi,xx0,yy0,xx,yy);
                     if (cosphi > 0 ) yy += ylabel;
                     if (cosphi < 0 ) yy -= ylabel;
                  }
               }
               if (!optionY || (x0 == x1)) {
                  if (!optionText) {
                     if (first > last)  strncpy(chtemp, " ", 256);
                     else               strncpy(chtemp, &label[first], 256);
                     if (fNModLabs) ChangeLabelAttributes(k+1, nlabels, textaxis, chtemp);
                     typolabel = chtemp;
                     if (!optionTime) typolabel.ReplaceAll("-", "#minus");
							//std::cout << "New hope1" << std::endl;
                     textaxis->PaintLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                           gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),
                           textaxis->GetTextAngle(),
                           textaxis->GetTextSize(),
                           typolabel.Data());

							//RADEK start
							//TLatex *lat = new TLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                           //gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()), typolabel.Data() );
							TLatex *lat = new TLatex(xx, yy, typolabel.Data());
							CopyLatexStyleNDC(lat, textaxis);
							lat->Draw();
							//std::cout << "Holcicka " << GetName() << std::endl;
							lat->SetName(TString::Format("LABEL%s", GetName()));
							recsNDC.push_back(GetNDCmy(lat));
							//if(!isVisual) lat->Delete();

                     if (fNModLabs) ResetLabelAttributes(textaxis);
                  }
                  else  {
                     if (optionText == 1) textaxis->PaintLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                                                   gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),
                                                   0,
                                                   textaxis->GetTextSize(),
                                                   fAxis->GetBinLabel(k+fAxis->GetFirst()));
                  }
               }
               else {

// Text alignment is down
                  if (!optionText)     lnlen = last-first+1;
                  else {
                     if (k+1 > nhilab) lnlen = 0;
                  }
                  for ( l=1; l<=lnlen; l++) {
                     if (!optionText) *chtemp = label[first+l-2];
                     else {
                        if (lnlen == 0) strncpy(chtemp, " ", 256);
                        else            strncpy(chtemp, "1", 256);
                     }
                     typolabel = chtemp;
                     typolabel.ReplaceAll("-", "#minus");
							std::cout << "New hope3" << std::endl;
                     textaxis->PaintLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                           gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),
                           0,
                           textaxis->GetTextSize(),
                           typolabel.Data());
                     yy -= charheight*1.3;
                  }
               }
            }

//   We use the format x 10 ** n

            if (flexe && !optionText && nexe)  {
               snprintf(label,256,"#times10^{%d}", nexe);
               if (x0 != x1) { xfactor = axis_length+0.1*charheight; yfactor = 0; }
               else          { xfactor = y1-y0+0.1*charheight; yfactor = 0; }
               Rotate (xfactor,yfactor,cosphi,sinphi,x0,y0,xx,yy);
               textaxis->SetTextAlign(11);
               if (GetLabelFont()%10 < 2) // force TLatex mode in PaintLatex
                  textaxis->SetTextFont((Int_t)(GetLabelFont()/10)*10+2);
               if (fAxis && !strcmp(fAxis->GetName(),"xaxis")) {
                  xx = xx + fXAxisExpXOffset;
                  yy = yy + fXAxisExpYOffset;
               }
               if (fAxis && !strcmp(fAxis->GetName(),"yaxis")) {
                  xx = xx + fYAxisExpXOffset;
                  yy = yy + fYAxisExpYOffset;
               }
               typolabel = label;
               typolabel.ReplaceAll("-", "#minus");
               textaxis->PaintLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                           gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),
                           0,
                           textaxis->GetTextSize(),
                           typolabel.Data());
					//std::cout << "PUSINKA " << typolabel << std::endl;
					TLatex *lat = new TLatex(xx, yy, typolabel.Data());
					lat->SetName("ExpInfo");
					CopyLatexStyleNDC(lat, textaxis);
					lat->Draw();
					//std::cout << "Holcicka " << GetName() << std::endl;
					//recsNDC.push_back(GetNDCmyold(lat, GetName()));


            }
         }
      }
   }

// Log axis

   if (optionLog && ndiv) {
      UInt_t xi1=0,xi2=0,wi=0,yi1=0,yi2=0,hi=0,xl=0,xh=0;
      Bool_t firstintlab = kTRUE, overlap = kFALSE;
      if ((wmin == wmax) || (ndiv == 0))  {
         Error(where, "wmin (%f) == wmax (%f), or ndiv == 0", wmin, wmax);
         goto L210;
      }
      if (wmin <= 0)   {
         Error(where, "negative logarithmic axis");
         goto L210;
      }
      if (wmax <= 0)     {
         Error(where, "negative logarithmic axis");
         goto L210;
      }
      xmnlog = TMath::Log10(wmin);
      if (xmnlog > 0) xmnlog += 1.E-6;
      else            xmnlog -= 1.E-6;
      x00    = 0;
      x11    = axis_length;
      h2     = TMath::Log10(wmax);
      h2sav  = h2;
      if (h2 > 0) h2 += 1.E-6;
      else        h2 -= 1.E-6;
      ih1    = int(xmnlog);
      ih2    = 1+int(h2);
      nbinin = ih2-ih1+1;
      axmul  = (x11-x00)/(h2sav-xmnlog);

// Plot decade and intermediate tick marks
      decade      = ih1-2;
      labelnumber = ih1;
      if ( xmnlog > 0 && (xmnlog-double_t(ih1) > 0) ) labelnumber++;
      for (j=1; j<=nbinin; j++) {

// Plot decade
         firstintlab = kTRUE, overlap = kFALSE;
         decade++;
         if (x0 == x1 && j == 1) ylabel += charheight*0.33;
         if (y0 == y1 && j == 1) ylabel -= charheight*0.65;
         xone = x00+axmul*(double_t(decade)-xmnlog);
         //the following statement is a trick to circumvent a gcc bug
         if (j < 0) printf("j=%d\n",j);
         if (x00 > xone) goto L160;
         if ((xone-x11)>epsilon) break;
         xtwo = xone;
         y    = 0;
         if (!mside) y -= atick[0];
         Rotate(xone,y,cosphi,sinphi,x0,y0,xpl2,ypl2);
         Rotate(xtwo,atick[0],cosphi,sinphi,x0,y0,xpl1,ypl1);
         if (optionVert) {
            if ((x0 != x1) && (y0 != y1)) {
               if (mside) {
                  xpl1=xpl2;
                  if (cosphi > 0) ypl1 = ypl2 + atick[0];
                  else            ypl1 = ypl2 - atick[0];
               }
               else {
                  xpl1 = 0.5*(xpl1 + xpl2);
                  xpl2 = xpl1;
                  ypl1 = 0.5*(ypl1 + ypl2) + atick[0];
                  ypl2 = 0.5*(ypl1 + ypl2) - atick[0];
               }
            }
         }
         if (!drawGridOnly) PaintLineNDC(xpl1, ypl1, xpl2, ypl2);

         if (optionGrid) {
            Rotate(xone,0,cosphi,sinphi,x0,y0,xpl2,ypl2);
            Rotate(xone,grid_side*gridlength,cosphi,sinphi,x0,y0,xpl1,ypl1);
            linegrid->PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
         }

         if (!drawGridOnly && !optionUnlab)  {

// We generate labels (numeric only).
            if (noExponent) {
               rlab = TMath::Power(10,labelnumber);
               snprintf(label,256, "%f", rlab);
               LabelsLimits(label,first,last);
               while (last > first) {
                  if (label[last] != '0') break;
                  label[last] = 0;
                  last--;
               }
               if (label[last] == '.') {label[last] = 0; last--;}
            } else {
               snprintf(label,256, "%d", labelnumber);
               LabelsLimits(label,first,last);
            }
            Rotate (xone,ylabel,cosphi,sinphi,x0,y0,xx,yy);
            if ((x0 == x1) && !optionPara) {
               if (lside < 0) {
                  if (mside < 0) {
                     if (labelnumber == 0) nch=1;
                     else                  nch=2;
                     xx    += nch*charheight;
                  } else {
                     xx += 0.25*charheight;
                  }
               }
               xx += 0.25*charheight;
            }
            if ((y0 == y1) && !optionDown && !optionUp) {
               if (noExponent) yy += 0.33*charheight;
            }
            if (n1a == 0)goto L210;
            kmod = nbinin/n1a;
            if (kmod == 0) kmod=1000000;
            if ((nbinin <= n1a) || (j == 1) || (j == nbinin) || ((nbinin > n1a)
            && (j%kmod == 0))) {

					std::cout << "RADEKhura " << std::endl;
               if (labelnumber == 0) {
                  textaxis->PaintTextNDC(xx,yy,"1");


						/*
						TText *tex = new TText(xx, yy, "1");
						CopyTextStyleNDC(tex, textaxis);
						tex->Draw();
						std::cout << "MMMMMMMMMMMM 1 : " << textaxis->GetTextAlign() << std::endl;
						unsigned w, h;
						tex->GetBoundingBox(w, h);
						std::cout << "AhojHolka " << w << " "<< h << std::endl;
						RectangleNDC Rec;
						Rec.fX = xx;
						Rec.fY = yy;
						Rec.fWidth = w / double(gPad->GetAbsWNDC()*gPad->GetWw());
						Rec.fHeight = h / double(gPad->GetAbsHNDC()*gPad->GetWh());
						if(!strcmp(GetName(), "xaxis")) {
							Rec.fX -= Rec.fWidth/2.;
						}
						else {
							Rec.fX -= Rec.fWidth*1.4;
							Rec.fY -= Rec.fHeight*0.5;

						}
						Rec.tex = tex;
						recsNDC.push_back(Rec);
						*/

						TLatex *lat = new TLatex(xx, yy, "1");
						CopyLatexStyleNDC(lat, textaxis);
						lat->Draw();
						lat->SetName(TString::Format("LABEL%s", GetName()));
						recsNDC.push_back(GetNDCmy(lat));



               } else if (labelnumber == 1) {
                  textaxis->PaintTextNDC(xx,yy,"10");

						/*
						TText *tex = new TText(xx, yy, "10");
						CopyTextStyleNDC(tex, textaxis);
						std::cout << "MMMMMMMMMMMM 10 : " << textaxis->GetTextAlign() << std::endl;
						tex->Draw();
						unsigned w, h;
						tex->GetBoundingBox(w, h);
						std::cout << "AhojHolka " << w << " "<< h << std::endl;
						RectangleNDC Rec;
						Rec.fX = xx;
						Rec.fY = yy;
						Rec.fWidth = w / double(gPad->GetAbsWNDC()*gPad->GetWw());
						Rec.fHeight = h / double(gPad->GetAbsHNDC()*gPad->GetWh());
						if(!strcmp(GetName(), "xaxis")) {
							Rec.fX -= Rec.fWidth/2.;
						}
						else {
							Rec.fX -= Rec.fWidth*1.0;
							Rec.fY -= Rec.fHeight*0.5;

						}
						Rec.tex = tex;
						recsNDC.push_back(Rec);
						*/

						TLatex *lat = new TLatex(xx, yy, "10");
						CopyLatexStyleNDC(lat, textaxis);
						lat->Draw();
						lat->SetName(TString::Format("LABEL%s", GetName()));
						recsNDC.push_back(GetNDCmy(lat));


               } else {
                  if (noExponent) {
                     textaxis->PaintTextNDC(xx,yy,&label[first]);
							std::cout << "RADEKhura nova3 PUSSSSAAAAAAAAAAAA " <<xx<<" "<<yy << std::endl;
							TLatex *lat = new TLatex(xx, yy, &label[first]);
							CopyLatexStyleNDC(lat, textaxis);
							lat->Draw();
							lat->SetName(TString::Format("LABEL%s", GetName()));
							recsNDC.push_back(GetNDCmy(lat));
                  } else {
                        snprintf(chtemp,256, "10^{%d}", labelnumber);
                        typolabel = chtemp;
                        typolabel.ReplaceAll("-", "#minus");
                        textaxis->PaintLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
                                             gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),
                                             0, textaxis->GetTextSize(), typolabel.Data());
								//cout << "RADEKhura nova4 " <<xx<<" "<<yy << std::endl;
								//TLatex *lat = static_cast<TLatex*>(textaxis->Clone());
								//lat->SetTextColor(kBlue);
								//lat->DrawLatexNDC(xx, yy, typolabel.Data());

								std::cout << "MMMMMMMMMMMM rest : " << textaxis->GetTextAlign() << std::endl;
								std::cout << "RADEKhura nova1 "<<xx<<" "<<yy << std::endl;
								TLatex *lat = new TLatex(xx, yy, typolabel.Data());
								CopyLatexStyleNDC(lat, textaxis);
								lat->Draw();
								lat->SetName(TString::Format("LABEL%s", GetName()));
								recsNDC.push_back(GetNDCmy(lat));
								//if(!isVisual) lat->Delete();

								//lat->DrawLatex(gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1()),
														  //gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1()),typolabel.Data());
                  }
               }
            }
            labelnumber++;
         }
L160:
         for (k=2;k<10;k++) {

// Plot intermediate tick marks
            xone = x00+axmul*(TMath::Log10(double_t(k))+double_t(decade)-xmnlog);
            if (x00 > xone) continue;
            if (xone > x11) goto L200;
            y = 0;
            if (!mside)  y -= atick[1];
            xtwo = xone;
            Rotate(xone,y,cosphi,sinphi,x0,y0,xpl2,ypl2);
            Rotate(xtwo,atick[1],cosphi,sinphi,x0,y0,xpl1,ypl1);
            if (optionVert) {
               if ((x0 != x1) && (y0 != y1)) {
                  if (mside) {
                     xpl1 = xpl2;
                     if (cosphi > 0) ypl1 = ypl2 + atick[1];
                     else            ypl1 = ypl2 - atick[1];
                  }
                  else {
                     xpl1 = 0.5*(xpl1+xpl2);
                     xpl2 = xpl1;
                     ypl1 = 0.5*(ypl1+ypl2) + atick[1];
                     ypl2 = 0.5*(ypl1+ypl2) - atick[1];
                  }
               }
            }
            idn = n1a*2;
            if ((nbinin <= idn) || ((nbinin > idn) && (k == 5))) {
               if (!drawGridOnly) PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
// Draw the intermediate LOG labels if requested

               if (moreLogLabels && !optionUnlab && !drawGridOnly && !overlap) {
                  if (noExponent) {
                     rlab = double_t(k)*TMath::Power(10,labelnumber-1);
                     snprintf(chtemp,256, "%g", rlab);
                  } else {
                     if (labelnumber-1 == 0) {
                        snprintf(chtemp,256, "%d", k);
                     } else if (labelnumber-1 == 1) {
                        snprintf(chtemp,256, "%d", 10*k);
                     } else {
                        snprintf(chtemp,256, "%d#times10^{%d}", k, labelnumber-1);
                     }
                  }
                  Rotate (xone,ylabel,cosphi,sinphi,x0,y0,xx,yy);
                  if ((x0 == x1) && !optionPara) {
                     if (lside < 0) {
                        if (mside < 0) {
                           if (labelnumber == 0) nch=1;
                           else                  nch=2;
                           xx    += nch*charheight;
                        } else {
                           if (labelnumber >= 0) xx    += 0.25*charheight;
                           else                  xx    += 0.50*charheight;
                        }
                     }
                     xx += 0.25*charheight;
                  }
                  if ((y0 == y1) && !optionDown && !optionUp) {
                     if (noExponent) yy += 0.33*charheight;
                  }
                  if (optionVert) {
                     if ((x0 != x1) && (y0 != y1)) {
                        Rotate(xone,ylabel,cosphi,sinphi,x0,y0,xx,yy);
                        if (cosphi > 0) yy += ylabel;
                        else            yy -= ylabel;
                     }
                  }
                  textaxis->SetTitle(chtemp);
                  double_t u = gPad->GetX1() + xx*(gPad->GetX2() - gPad->GetX1());
                  double_t v = gPad->GetY1() + yy*(gPad->GetY2() - gPad->GetY1());
                  if (firstintlab) {
                     textaxis->GetBoundingBox(wi, hi); wi=(UInt_t)(wi*1.3); hi=(UInt_t)(hi*1.3);
                     xi1 = gPad->XtoAbsPixel(u);
                     yi1 = gPad->YtoAbsPixel(v);
                     firstintlab = kFALSE;
                     typolabel = chtemp;
                     typolabel.ReplaceAll("-", "#minus");
                     textaxis->PaintLatex(u,v,0,textaxis->GetTextSize(),typolabel.Data());

							//RADEK
							TLatex *lat =(TLatex*) textaxis->Clone();
							lat->SetTextColor(textaxis->GetTextColor());
							lat->SetNDC();
							lat->SetText(xx, yy, typolabel.Data());
							lat->Draw();
							lat->SetName(TString::Format("LABEL%s", GetName()));
							recsNDC.push_back(GetNDCmy(lat));
							//if(!isVisual) lat->Delete();


                  } else {
                     xi2 = gPad->XtoAbsPixel(u);
                     yi2 = gPad->YtoAbsPixel(v);
                     xl = TMath::Min(xi1,xi2);
                     xh = TMath::Max(xi1,xi2);
                     if ((x0 == x1 && yi1-hi <= yi2) || (y0 == y1 && xl+wi >= xh)){
                        overlap = kTRUE;
                     } else {
                        xi1 = xi2;
                        yi1 = yi2;
                        textaxis->GetBoundingBox(wi, hi); wi=(UInt_t)(wi*1.3); hi=(UInt_t)(hi*1.3);
                        typolabel = chtemp;
                        typolabel.ReplaceAll("-", "#minus");
                        textaxis->PaintLatex(u,v,0,textaxis->GetTextSize(),typolabel.Data());
								//RADEK
								TLatex *lat =(TLatex*) textaxis->Clone();
								lat->SetTextColor(textaxis->GetTextColor());
								lat->SetNDC();
								lat->SetText(xx, yy, typolabel.Data());
								lat->Draw();
								lat->SetName(TString::Format("LABEL%s", GetName()));
								recsNDC.push_back(GetNDCmy(lat));
								//if(!isVisual) lat->Delete();
								//std::cout <<"TEREZKA2 "<<u<<" "<<v << std::endl;
                     }
                  }
               }

// Draw the intermediate LOG grid if only three decades are requested
               if (optionGrid && nbinin <= 5 && ndiv > 100) {
                  Rotate(xone,0,cosphi,sinphi,x0,y0,xpl2, ypl2);
                  Rotate(xone,grid_side*gridlength,cosphi,sinphi,x0,y0, xpl1,ypl1);
                  linegrid->PaintLineNDC(xpl1, ypl1, xpl2, ypl2);
               }
            }  //endif ((nbinin <= idn) ||
         }  //endfor (k=2;k<10;k++)
      } //endfor (j=1; j<=nbinin; j++)
L200:
      Int_t dummy = 0; if (dummy) { }
   }  //endif (optionLog && ndiv)


L210:
   if (optionGrid) delete linegrid;
   delete textaxis;
}



};



inline double RemoveOverlaps(TVirtualPad *pad, TAxis *ax, bool remFirst, bool remLast, bool isVisual = false)
{
	//std::cout << "Helenka " << ax->GetName() << std::endl;
	pad->cd();
	pad->Update();

	bool isX;
	if(!strcmp(ax->GetName(), "xaxis"))
		isX = true;
	else if(!strcmp(ax->GetName(), "yaxis"))
		isX = false;
	else
		return 0;

	double valMin, valMax;
	//int PxMin, PxMax;

	if(isX) {
		valMin = pad->GetUxmin();
		valMax = pad->GetUxmax();
		//PxMin = pad->XtoPixel(valMin);
		//PxMax = pad->XtoPixel(valMax);
	}
	else {
		valMin = pad->GetUymin();
		valMax = pad->GetUymax();
		//PxMin = pad->YtoPixel(valMin);
		//PxMax = pad->YtoPixel(valMax);
	}

	if((isX && pad->GetLogx()) || (!isX && pad->GetLogy()) ) {
		valMin = pow(10, valMin);
		valMax = pow(10, valMax);

	}

	myTGaxis *gAx = new myTGaxis();
	gAx->SetName(ax->GetName());
	gAx->isVisual = isVisual;
	gAx->SetLabelFont(ax->GetLabelFont());
	gAx->SetLabelSize(ax->GetLabelSize());
	gAx->SetLabelOffset(ax->GetLabelOffset());
	gAx->SetTitle(ax->GetTitle());
	gAx->SetTitleOffset(ax->GetTitleOffset());
	gAx->SetTitleSize(ax->GetTitleSize());
	gAx->SetTitleFont(ax->GetTitleFont());

	gAx->SetLineColor(kRed);

	ax->SetLabelSize(0);


	int ndiv = ax->GetNdivisions();
	//double wMin = gPad->GetUxmin(), wMax = gPad->GetUxmax();

	gAx->SetMoreLogLabels(ax->GetMoreLogLabels());
	gAx->SetNoExponent(ax->GetNoExponent());
	gAx->SetDecimals(ax->GetDecimals());
	gAx->CenterTitle(ax->GetCenterTitle());

	if(isX) {
		if(!pad->GetLogx())
			gAx->PaintAxis(pad->GetUxmin(), pad->GetUymin(), pad->GetUxmax(), pad->GetUymin(), valMin, valMax, ndiv, "");
		else {
			std::cout << "RADEKOUTPUT " << pow(10,pad->GetUxmin()) <<" "<< pow(10,pad->GetUxmax()) << std::endl;
			std::cout << "RADEKOUTPUT " << valMin <<" "<< valMax << std::endl;
			//gPad->SetLogx();
			//gAx->PaintAxis(pow(10,pad->GetUxmin()), pad->GetUymin(), pow(10,pad->GetUxmax()), pad->GetUymin(), valMin, valMax, ndiv, "G");
			gAx->PaintAxis(pad->GetUxmin(), pad->GetUymin(), pad->GetUxmax(), pad->GetUymin(), valMin, valMax, ndiv, "G");
			//gAx->SetLabelColor(kRed);
			//gAx->DrawAxis(pow(10,pad->GetUxmin()), pad->GetUymin(), pow(10,pad->GetUxmax()), pad->GetUymin(), valMin, valMax, ndiv, "G");
		}
	}
	else {
		if(!pad->GetLogy())
			gAx->PaintAxis(pad->GetUxmin(), pad->GetUymin(), pad->GetUxmin(), pad->GetUymax(), valMin, valMax, ndiv, "");
		else
			gAx->PaintAxis(pad->GetUxmin(), pad->GetUymin(), pad->GetUxmin(), pad->GetUymax(), valMin, valMax, ndiv, "G");
	}



	TList *li =  gPad->GetListOfPrimitives();
	for(unsigned i = 0; i < gAx->recsNDC.size(); ++i) {
			TLine *lin = new TLine;
			lin->SetLineColor(kRed);
			double myX1 = gAx->recsNDC[i].fX;
			double myX2 = gAx->recsNDC[i].fX+ gAx->recsNDC[i].fWidth;
			double myY1 = gAx->recsNDC[i].fY;
			double myY2 = gAx->recsNDC[i].fY+ gAx->recsNDC[i].fHeight;
			if(isVisual) {
				lin->DrawLineNDC(myX1, myY1, myX2, myY1);
				lin->DrawLineNDC(myX1, myY2, myX2, myY2);
				lin->DrawLineNDC(myX1, myY1, myX1, myY2);
				lin->DrawLineNDC(myX2, myY1, myX2, myY2);
			}
			//std::cout << "Marketka "<<isX<<" "<<i<<" "<<   myX1 <<" "<< myX2 <<" : "<< myY1<<" "<<myY2<< std::endl;

			bool isOutside = (isX && ((remFirst && myX1 < 0) || (remLast  && myX2 > 1))) ||
			                (!isX && ((remFirst && myY1 < 0) || (remLast && myY2 > 1) ));

			if(isOutside) {
				//ax->ChangeLabel(i+1, -1, 0);
				if(gAx->recsNDC[i].lat) li->Remove(gAx->recsNDC[i].lat);
				if(gAx->recsNDC[i].tex) li->Remove(gAx->recsNDC[i].tex);
				//std::cout << "Removing " << i+1 << std::endl;
			}


	}




	//Plot title rect
	TLine *lin = new TLine;
	lin->SetLineColor(kBlue);

	double myX1 = gAx->rTitleNDC.fX;
	double myX2 = gAx->rTitleNDC.fX+ gAx->rTitleNDC.fWidth;
	double myY1 = gAx->rTitleNDC.fY;
	double myY2 = gAx->rTitleNDC.fY+ gAx->rTitleNDC.fHeight;
	if(isVisual) {
		lin->DrawLineNDC(myX1, myY1, myX2, myY1);
		lin->DrawLineNDC(myX1, myY2, myX2, myY2);
		lin->DrawLineNDC(myX1, myY1, myX1, myY2);
		lin->DrawLineNDC(myX2, myY1, myX2, myY2);
	}

	/*
	myX1 = gAx->rTitleMyNDC.fX;
	myX2 = gAx->rTitleMyNDC.fX+ gAx->rTitleMyNDC.fWidth;
	myY1 = gAx->rTitleMyNDC.fY;
	myY2 = gAx->rTitleMyNDC.fY+ gAx->rTitleMyNDC.fHeight;
	if(isVisual) {
		lin->DrawLineNDC(myX1, myY1, myX2, myY1);
		lin->DrawLineNDC(myX1, myY2, myX2, myY2);
		lin->DrawLineNDC(myX1, myY1, myX1, myY2);
		lin->DrawLineNDC(myX2, myY1, myX2, myY2);
	}
	*/



	gPad->Update();


	if(strlen(gAx->GetTitle()) > 0)
		return gAx->GetMinimalOffset();
	else
		return 0;


}

inline void RemoveLabels(TVirtualPad *pad, const char *axis)
{
	//pad->GetListOfPrimitives()->Print();
	TList *li =  pad->GetListOfPrimitives();
	for( auto && p : *li) {
		TString name = p->GetName();
		//std::cout << "RADECEK " << name << std::endl;
		if(name == TString("LABEL")+axis) {
			//std::cout << "Removing " << name <<" "<< p->ClassName()<< std::endl;
			li->Remove(p);
		}
		//if(gAx->recsNDC[i].tex) li->Remove(gAx->recsNDC[i].tex);
	}
	pad->Modified();
}



inline int GetIntersect(Edge e, Point p, double dirX, double dirY, double &shift1, double &shift2)
{
	const double eps = 1e-9;
  double a11 = e.p2.x - e.p1.x; 
  double a21 = e.p2.y - e.p1.y; 
	
  double a12 = -dirX; 
  double a22 = -dirY; 

  double b1 = p.x - e.p1.x;
  double b2 = p.y - e.p1.y;


  double det = a11*a22 - a12*a21;

  double detE = b1*a22 - a12*b2; //Edge
  double detS = a11*b2 - b1*a21; //Shift

  if( abs(det) > eps ) {
    double ePos = detE/det;
    if(ePos < 0 || ePos > 1)
      return 0;
    shift1 = detS / det;
    return 1;
  }
  else {
    if( abs(detE) > eps || abs(detS) > eps )
      return 0;
    if( abs(dirX) > abs(dirY) ) {
      shift1 = (e.p1.x - p.x ) / dirX;
      shift2 = (e.p2.x - p.x ) / dirX;
		assert(e.p1.x == e.p1.x);
		assert(p.x == p.x);
		assert(dirX == dirX);
      return 2;
    }
    else {
      shift1 = (e.p1.y - p.y ) / dirY;
      shift2 = (e.p2.y - p.y ) / dirY;
		assert(e.p1.y == e.p1.y);
		assert(p.y == p.y);
		//std::cout <<"Marketka "<< dirX << " "<< dirY << std::endl;
		assert(dirY == dirY);
      return 2;
    }
  }

}



double myTGaxis::GetMinimalOffset()
{
	
	//vector<RectangleNDC> recsNDC;
	//RectangleNDC rTitleNDC;
	//RectangleNDC rTitleMyNDC;

	auto AddRectangle = [](const RectangleNDC &rec, std::vector<Point> &Points, std::vector<Edge> &Edges) {
		Point P1(rec.fX, rec.fY);
		Point P2(rec.fX+rec.fWidth, rec.fY);
		Point P3(rec.fX, rec.fY+rec.fHeight);
		Point P4(rec.fX+rec.fWidth, rec.fY+rec.fHeight);

		Points.push_back(P1);
		Points.push_back(P2);
		Points.push_back(P3);
		Points.push_back(P4);

		Edge e1, e2, e3, e4;
		e1.p1 = P1; e1.p2 = P2;
		e2.p1 = P2; e2.p2 = P3;
		e3.p1 = P3; e3.p2 = P4;
		e4.p1 = P4; e3.p2 = P1;
		Edges.push_back(e1);
		Edges.push_back(e2);
		Edges.push_back(e3);
		Edges.push_back(e4);
	};

	std::vector<Point> tPoints, tPointsMy;
	std::vector<Edge>  tEdges, tEdgesMy;
	AddRectangle(rTitleNDC, tPoints, tEdges);
	AddRectangle(rTitleMyNDC, tPointsMy, tEdgesMy);
	double dirX = tPointsMy[0].x - tPoints[0].x;
	double dirY = tPointsMy[0].y - tPoints[0].y;
	double l = hypot(dirX, dirY);
	dirX /= l;
	dirY /= l;

	std::vector<Edge> lEdges;
	std::vector<Point> lPoints;
	for(const RectangleNDC &r : recsNDC) {
		AddRectangle(r, lPoints, lEdges);
	}

	double MinL = 1e30, MinT = -1e30, MinT0 = -1e30;
	for(const Point &p : lPoints) {
		MinL = abs(dirX)>abs(dirY) ? std::min(MinL, p.x) : std::min(MinL, p.y);
		//std::cout << "Label " << p.y << std::endl;
	}
	for(const Point &p : tPoints) {
		MinT = abs(dirX)>abs(dirY) ? std::max(MinT, p.x) : std::max(MinT, p.y);
		//std::cout << "PointOrg " << p.y << std::endl;
	}
	for(const Point &p : tPointsMy) {
		MinT0 = abs(dirX)>abs(dirY) ? std::max(MinT0, p.x) : std::max(MinT0, p.y);
		//std::cout << "PointMy " << p.y << std::endl;
	}

	double Factor = (MinL - MinT0) / (MinT - MinT0);
	//std::cout << "Holkaaa " << MinT <<" "<< MinT0 << std::endl;
	//std::cout << "MinL " << MinL  << std::endl;

	/*
	vector<double> intersects;
	for(const Point &p : tPoints) {
		double shift1, shift2;
		for(Edge &e : lEdges) {
			std::cout << "Kamila " << dirX <<" "<< dirY << std::endl;
			int status = GetIntersect(e, p, dirX, dirY, shift1, shift2);
			if(status == 2) {
				intersects.push_back(shift1);
				intersects.push_back(shift2);
				assert(shift1 == shift1);
				assert(shift2 == shift2);
			}
			else if(status == 1) {
				intersects.push_back(shift1);
			}
		}
	}
	for(const Point &p : lPoints) {
		double shift1, shift2;
		for(Edge &e : tEdges) {
			int status = GetIntersect(e, p, -dirX, -dirY, shift1, shift2);
			if(status == 2) {
				intersects.push_back(shift1);
				intersects.push_back(shift2);
			}
			else if(status == 1) {
				intersects.push_back(shift1);
			}
		}
	}
	double maxInter = *max_element(intersects.begin(), intersects.end());

	double myOffsetEstimate = (1 + maxInter/l * 0.2) * GetTitleOffset();
	*/
	double myOffsetEstimate = abs(Factor) * GetTitleOffset();

	//std::cout << "directions " << dirX << " "<< dirY << std::endl;
	//std::cout << "myNumber " << myOffsetEstimate << std::endl;
	return myOffsetEstimate;

}
