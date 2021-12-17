//
// Split View
//

#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TColor.h>
#include <TRandom.h>
#include <TMath.h>
#include <TEveBox.h>
#include <TEveArrow.h>
#include <TEveBrowser.h>

class EvNavHandler
{
public:
  void Fwd()
  {

      printf("Already at last event.\n");

  }
  void Bck()
  {

      printf("Already at first event.\n");

  }
};

void make_gui()
{
  // Create minimal GUI for event navigation.
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);
  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
  {
    TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    TGPictureButton* b = 0;
    EvNavHandler    *fh = new EvNavHandler;
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
    hf->AddFrame(b);
    b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
    hf->AddFrame(b);
    b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
  }
  frmMain->AddFrame(hf);
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);
}


void run_vis_3()
{
  TEveManager::Create();

  gEve->GetBrowser()->GetTabRight()->SetTab(1);
  make_gui();
  
  auto ps = new TEvePointSet();
  
  ps->SetNextPoint(1.5, 1.5, 0.1);
  ps->SetNextPoint(7.6, 1.3, 0.1);
  ps->SetNextPoint(1.7, 7.1, 0.1);
  ps->SetNextPoint(7.8, 6.9, 0.1);

  
  ps->SetMarkerColor(3);
  ps->SetMarkerSize(2);
  ps->SetMarkerStyle(4);
  
  double origin[3] = {0, 0, 0};
  
  TEveBox* b0 = new TEveBox();
  b0->SetVertex(0, 0, 0, 0);
  b0->SetVertex(1, 4, 0, 0);
  b0->SetVertex(2, 4, 4, 0);
  b0->SetVertex(3, 0, 4, 0);
  b0->SetVertex(4, 0, 0, 0.03);
  b0->SetVertex(5, 4, 0, 0.03);
  b0->SetVertex(6, 4, 4, 0.03);
  b0->SetVertex(7, 0, 4, 0.03);
  
  TEveBox* b1 = new TEveBox();
  b1->SetVertex(0, 6, 0, 0);
  b1->SetVertex(1, 10, 0, 0);
  b1->SetVertex(2, 10, 4, 0);
  b1->SetVertex(3, 6, 4, 0);
  b1->SetVertex(4, 6, 0, 0.03);
  b1->SetVertex(5, 10, 0, 0.03);
  b1->SetVertex(6, 10, 4, 0.03);
  b1->SetVertex(7, 6, 4, 0.03);
  
  TEveBox* b2 = new TEveBox();
  b2->SetVertex(0, 0, 6, 0);
  b2->SetVertex(1, 4, 6, 0);
  b2->SetVertex(2, 4, 10, 0);
  b2->SetVertex(3, 0, 10, 0);
  b2->SetVertex(4, 0, 6, 0.03);
  b2->SetVertex(5, 4, 6, 0.03);
  b2->SetVertex(6, 4, 10, 0.03);
  b2->SetVertex(7, 0, 10, 0.03);
  
  TEveBox* b3 = new TEveBox();
  b3->SetVertex(0, 6, 6, 0);
  b3->SetVertex(1, 10, 6, 0);
  b3->SetVertex(2, 10, 10, 0);
  b3->SetVertex(3, 6, 10, 0);
  b3->SetVertex(4, 6, 6, 0.03);
  b3->SetVertex(5, 10, 6, 0.03);
  b3->SetVertex(6, 10, 10, 0.03);
  b3->SetVertex(7, 6, 10, 0.03);
  
  gEve->AddElement(b0);
  gEve->AddElement(b1);
  gEve->AddElement(b2);
  gEve->AddElement(b3);

  
  gEve->AddElement(ps);
  
  gEve->Redraw3D();
}
