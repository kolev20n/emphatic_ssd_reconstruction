#include <TGFrame.h>
#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>

class MyMainFrame : public TGMainFrame
{
private:
  TGMainFrame *fMain;
  TRootEmbeddedCanvas *fEcanvas;
public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame()
  {
    fMain->Cleanup();
    delete fMain;
  }
  void DoDraw()
  {
    TF1 *f1 = new TF1("f1","sin(x)/x",0,gRandom->Rndm()*10); f1->SetLineWidth(3);
    f1->Draw();
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
  }
  //ClassDef(MyMainFrame,0)
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h)
: TGMainFrame(p,w,h)
{
  // Creates widgets of the example
  fEcanvas = new TRootEmbeddedCanvas ("Ecanvas",this,200,200);
  AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  TGHorizontalFrame *hframe=new TGHorizontalFrame(this, 200,40);
  TGTextButton *draw = new TGTextButton(hframe,"&Draw");
  draw->Connect("Clicked()","MyMainFrame",this,"DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
  TGTextButton *exit = new TGTextButton(hframe,"&Exit ","gApplication->Terminate()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
  AddFrame(hframe,new TGLayoutHints(kLHintsCenterX,2,2,2,2));                                                                                                 SetWindowName("Simple Example");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

void example()
{
  new MyMainFrame(gClient->GetRoot(),200,200);
}

int main(int argc, char **argv)
{
  TApplication theApp("App",&argc,argv);
  example();
  theApp.Run();
  return 0;
}
