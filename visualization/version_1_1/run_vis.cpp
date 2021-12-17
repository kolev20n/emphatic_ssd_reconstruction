#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TColor.h>
#include <TRandom.h>
#include <TMath.h>
#include <TEveBox.h>
#include <TEveArrow.h>
#include <TEveBrowser.h>
#include <TGButton.h>
#include <TGTab.h>
#include <TSystem.h>
#include <iostream>
#include <TRint.h>

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
  frmMain->SetWindowName("YY GUI");
  frmMain->SetCleanup(kDeepCleanup);
  
  std::cout << "111" << std::endl;
  
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


int main(int argc, char **argv)
{
  std::cout << "000" << std::endl;
  
  TRint *app = new TRint("App", &argc, argv);
  
  TEveManager::Create();
  app->Run(kTRUE);
  
  gEve->GetBrowser()->GetTabRight()->SetTab(1);
  
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);
  
  std::cout << "111" << std::endl;
  
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
  
  auto ps = new TEvePointSet();
  
  ps->SetNextPoint(25, 30, 48);
  ps->SetNextPoint(25, 80, 50);
  ps->SetNextPoint(50, 30, 80);
  ps->SetNextPoint(90, 30, 48);
  ps->SetNextPoint(90, 10, 20);

  ps->SetMarkerColor(3);
  ps->SetMarkerSize(2);
  ps->SetMarkerStyle(4);
  
  double origin[3] = {0, 0, 0};
  
  TEveBox* b = new TEveBox();
  b->SetVertex(0, 0, 0, 0);
  b->SetVertex(1, 0, 40, 0);
  b->SetVertex(2, 0, 40, 40);
  b->SetVertex(3, 0, 0, 40);
  b->SetVertex(4, 2, 0, 0);
  b->SetVertex(5, 2, 40, 0);
  b->SetVertex(6, 2, 40, 40);
  b->SetVertex(7, 2, 0, 40);
  
  
  
  
  gEve->AddElement(b);
  
  gEve->AddElement(ps);
  
  
  
  gEve->Redraw3D();
  
  //app->Terminate(0);
  
  return 0;
}
