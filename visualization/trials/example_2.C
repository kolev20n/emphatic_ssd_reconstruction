#include <fstream>
#include <string>
#include <iostream>

#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TMath.h>
#include <Math/Vector3D.h>

using namespace std;
using namespace ROOT::Math;

class particle
{
public:
  int id;
  int code;
  double ekin_vertex;
  int n_plates;
  int plate_index[13];
  XYZVector plate_hits[13];
};

class cluster
{
public:
  int plate_index;
  int n_clusters;
  XYZVector points[100];
  XYZVector directions[100];
};

class track_option
{
public:
  int index_6, index_7, index_8;
  int strip;
  double distance;
  XYZVector point_on_d_plate;
  XYZVector track_line_point;
  XYZVector track_line_direction;
};

class MyMainFrame {
  RQ_OBJECT("MyMainFrame")
private:
  TGMainFrame         *fMain;
  TRootEmbeddedCanvas *ec1, *ec2, *ec3, *ec4, *ec5, *ec6, *ec7, *ec8, *ec9;
  
  ifstream in_file;
  
public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();
  void DoDraw();
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h)
{
  //in_file.open("temp_debug_midstream_downstream_with_tracks_with_field_all.txt", ios::in);
  //in_file.open("temp_debug_midstream_downstream_with_tracks_all.txt", ios::in);
  in_file.open("temp_debug_output.txt", ios::in);

  // Create a main frame
  fMain = new TGMainFrame(p,w,h);
  
  // Create canvas widget

  TGVerticalFrame *f3 = new TGVerticalFrame(fMain,1200,2000);
  TGHorizontalFrame *f1 = new TGHorizontalFrame(f3,1200,800);
  TGHorizontalFrame *f2 = new TGHorizontalFrame(f3,1200,800);

  ec1 = new TRootEmbeddedCanvas("Ecanvas1",f1,1200,600);
  ec2 = new TRootEmbeddedCanvas("Ecanvas2",f1,1200,600);
  ec3 = new TRootEmbeddedCanvas("Ecanvas3",f1,1200,600);
  ec4 = new TRootEmbeddedCanvas("Ecanvas4",f1,1200,600);
  ec5 = new TRootEmbeddedCanvas("Ecanvas5",f1,1200,600);
  ec6 = new TRootEmbeddedCanvas("Ecanvas6",f2,1200,600);
  ec7 = new TRootEmbeddedCanvas("Ecanvas7",f2,1200,600);
  ec8 = new TRootEmbeddedCanvas("Ecanvas8",f2,1200,600);
  ec9 = new TRootEmbeddedCanvas("Ecanvas9",f2,1200,600);

  f1->AddFrame(ec1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f1->AddFrame(ec2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f1->AddFrame(ec3, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f1->AddFrame(ec4, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f1->AddFrame(ec5, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f2->AddFrame(ec6, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f2->AddFrame(ec7, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f2->AddFrame(ec8, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f2->AddFrame(ec9, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));

  f3->AddFrame(f1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f3->AddFrame(f2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));

  fMain->AddFrame(f3, new TGLayoutHints(kLHintsExpandX |
                                              kLHintsExpandY, 10,10,10,1));

  // Create a horizontal frame widget with buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
  TGTextButton *draw = new TGTextButton(hframe,"&Draw");
  draw->Connect("Clicked()","MyMainFrame",this,"DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX,
                                           5,5,3,4));
  TGTextButton *exit = new TGTextButton(hframe,"&Exit",
                                        "gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,
                                           5,5,3,4));
  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
                                            2,2,2,2));
  
  // Set a name to the main frame
  fMain->SetWindowName("Simple Example");
  
  // Map all subwindows of main frame
  fMain->MapSubwindows();
  
  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());
  
  // Map main frame
  fMain->MapWindow();
}
void MyMainFrame::DoDraw()
{
  string card;
  int event_class;
  int cluster_multiplicity[13];
  double x, y, z;
  XYZVector mc_vertex, guessed_vertex;
  int n_particles = 0;
  particle the_particles[100];
  int n_cluster_plates = 0;
  cluster the_clusters[100];
  int n_track_options;
  track_option the_track_options[1000];
  int n_lines[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  TLine lines[100][9];
  const double sqrt2 = 0.707106781186548;
  double t;
  int n_tracks;
  int n_track_hits;
  int plate_index;
  
  // Draws function graphics in randomly chosen interval
  TF1 *f1 = new TF1("f1","sin(x)/x",0,gRandom->Rndm()*10);
  f1->SetLineWidth(3);
  
  TH2D *h41 = new TH2D("particles_plate_4", "particles_plate_4", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h51 = new TH2D("particles_plate_5", "particles_plate_5", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h61 = new TH2D("particles_plate_6", "particles_plate_6", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h71 = new TH2D("particles_plate_7", "particles_plate_7", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h81 = new TH2D("particles_plate_8", "particles_plate_8", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h91 = new TH2D("particles_plate_9", "particles_plate_9", 640, -3.84, 3.84, 640, -3.84, 3.84);
  TH2D *h101 = new TH2D("particles_plate_10", "particles_plate_10", 640, -3.84, 3.84, 640, -3.84, 3.84);
  TH2D *h111 = new TH2D("particles_plate_11", "particles_plate_11", 640, -3.84, 3.84, 640, -3.84, 3.84);
  TH2D *h121 = new TH2D("particles_plate_12", "particles_plate_12", 640, -3.84, 3.84, 640, -3.84, 3.84);

  TH2D *h62 = new TH2D("track_options", "track_options", 640, -1.92, 1.92, 640, -1.92, 1.92);
  
  TH2D *h93 = new TH2D("tracks", "tracks", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h103 = new TH2D("tracks", "tracks", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h113 = new TH2D("tracks", "tracks", 640, -1.92, 1.92, 640, -1.92, 1.92);
  TH2D *h123 = new TH2D("tracks", "tracks", 640, -1.92, 1.92, 640, -1.92, 1.92);

  h41->SetMarkerColor(4);
  h41->SetMarkerStyle(20);
  h41->SetMarkerSize(1.5);

  h51->SetMarkerColor(4);
  h51->SetMarkerStyle(20);
  h51->SetMarkerSize(1.5);
  
  h61->SetMarkerColor(4);
  h61->SetMarkerStyle(20);
  h61->SetMarkerSize(1.5);
  
  h71->SetMarkerColor(4);
  h71->SetMarkerStyle(20);
  h71->SetMarkerSize(1.5);
  
  h81->SetMarkerColor(4);
  h81->SetMarkerStyle(20);
  h81->SetMarkerSize(1.5);

  h91->SetMarkerColor(4);
  h91->SetMarkerStyle(20);
  h91->SetMarkerSize(1.5);
  
  h101->SetMarkerColor(4);
  h101->SetMarkerStyle(20);
  h101->SetMarkerSize(1.5);
  
  h111->SetMarkerColor(4);
  h111->SetMarkerStyle(20);
  h111->SetMarkerSize(1.5);
  
  h121->SetMarkerColor(4);
  h121->SetMarkerStyle(20);
  h121->SetMarkerSize(1.5);
  
  h62->SetMarkerColor(2);
  h62->SetMarkerStyle(43);
  h62->SetMarkerSize(2.0);
  
  h93->SetMarkerColor(3);
  h93->SetMarkerStyle(47);
  h93->SetMarkerSize(1.5);
  
  h103->SetMarkerColor(3);
  h103->SetMarkerStyle(47);
  h103->SetMarkerSize(1.5);
  
  h113->SetMarkerColor(3);
  h113->SetMarkerStyle(47);
  h113->SetMarkerSize(1.5);
  
  h123->SetMarkerColor(3);
  h123->SetMarkerStyle(47);
  h123->SetMarkerSize(1.5);
  
  TCanvas *fCanvas1 = ec1->GetCanvas();
  TCanvas *fCanvas2 = ec2->GetCanvas();
  TCanvas *fCanvas3 = ec3->GetCanvas();
  TCanvas *fCanvas4 = ec4->GetCanvas();
  TCanvas *fCanvas5 = ec5->GetCanvas();
  TCanvas *fCanvas6 = ec6->GetCanvas();
  TCanvas *fCanvas7 = ec7->GetCanvas();
  TCanvas *fCanvas8 = ec8->GetCanvas();
  TCanvas *fCanvas9 = ec9->GetCanvas();

  while (in_file >> card)
  {
    if (card == "event")
    {
      in_file >> event_class;
      for (int i = 0; i < 13; i++)
      {
        in_file >> cluster_multiplicity[i];
      }
      in_file >> x >> y >> z;
      mc_vertex.SetCoordinates(x, y, z);
      in_file >> x >> y >> z;
      guessed_vertex.SetCoordinates(x, y, z);
    }
    else if (card == "particle")
    {
      in_file >> the_particles[n_particles].id;
      in_file >> the_particles[n_particles].code;
      in_file >> the_particles[n_particles].ekin_vertex;
      in_file >> the_particles[n_particles].n_plates;
      for (int i = 0; i < the_particles[n_particles].n_plates; i++)
      {
        in_file >> the_particles[n_particles].plate_index[i];
        in_file >> x >> y >> z;
        the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
      }
      n_particles++;
    }
    else if (card == "clusters")
    {
      in_file >> the_clusters[n_cluster_plates].plate_index;
      in_file >> the_clusters[n_cluster_plates].n_clusters;
      for (int i = 0; i < the_clusters[n_cluster_plates].n_clusters; i++)
      {
        in_file >> x >> y >> z;
        the_clusters[n_cluster_plates].points[i].SetCoordinates(x, y, z);
        in_file >> x >> y >> z;
        the_clusters[n_cluster_plates].directions[i].SetCoordinates(x, y, z);
      }
      n_cluster_plates++;
    }
    else if (card == "track_options")
    {
      in_file >> n_track_options;
      for (int i = 0; i < n_track_options; i++)
      {
        in_file >> the_track_options[i].index_7 >> the_track_options[i].index_8 >> the_track_options[i].index_6 >> the_track_options[i].strip >> the_track_options[i].distance;
        in_file >> x >> y >> z;
        the_track_options[i].point_on_d_plate.SetCoordinates(x, y, z);
        in_file >> x >> y >> z;
        the_track_options[i].track_line_point.SetCoordinates(x, y, z);
        in_file >> x >> y >> z;
        the_track_options[i].track_line_direction.SetCoordinates(x, y, z);
      }
    }
    else if (card == "tracks")
    {
      in_file >> n_tracks;
      for (int i = 0; i < n_tracks; i++)
      {
        in_file >> n_track_hits;
        for (int j = 0; j < n_track_hits; j++)
        {
          in_file >> plate_index >> x >> y >> z;
          if (plate_index == 9)
          {
            h93->Fill(x, y);
          }
          else if (plate_index == 10)
          {
            h103->Fill(x, y);
          }
          else if (plate_index == 11)
          {
            h113->Fill(x, y);
          }
          else if (plate_index == 12)
          {
            h123->Fill(x, y);
          }
        }
      }
    }
    else if (card == "endevent")
    {
      for (int i = 0; i < n_particles; i++)
      {
        for (int j = 0; j < the_particles[i].n_plates; j++)
        {
          if (the_particles[i].plate_index[j] == 4)
          {
            h41->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 5)
          {
            h51->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 6)
          {
            h61->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 7)
          {
            h71->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 8)
          {
            h81->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 9)
          {
            h91->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 10)
          {
            h101->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 11)
          {
            h111->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
          else if (the_particles[i].plate_index[j] == 12)
          {
            h121->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
          }
        }
      }
      
      for (int i = 0; i < n_cluster_plates; i++)
      {
        for (int j = 0; j < the_clusters[i].n_clusters; j++)
        {
          if (the_clusters[i].plate_index == 4)
          {
            lines[j][i].SetX1(-1.92);
            lines[j][i].SetX2(1.92);
            lines[j][i].SetY1(the_clusters[i].points[j].Y());
            lines[j][i].SetY2(the_clusters[i].points[j].Y());
            n_lines[0]++;
          }
          if (the_clusters[i].plate_index == 5)
          {
            lines[j][i].SetY1(-1.92);
            lines[j][i].SetY2(1.92);
            lines[j][i].SetX1(the_clusters[i].points[j].X());
            lines[j][i].SetX2(the_clusters[i].points[j].X());
            n_lines[1]++;
          }
          if (the_clusters[i].plate_index == 6)
          {
            if (the_clusters[i].points[j].X() < 0.)
            {
              t = (-1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
            
              lines[j][i].SetX1(-1.92);
              lines[j][i].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
            
              t = (-1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i].SetY2(-1.92);
              n_lines[2]++;
            }
            else
            {
              t = (1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i].SetX1(1.92);
              lines[j][i].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i].SetY2(1.92);
              n_lines[2]++;
            }
          }
          if (the_clusters[i].plate_index == 7)
          {
            lines[j][i].SetX1(-1.92);
            lines[j][i].SetX2(1.92);
            lines[j][i].SetY1(the_clusters[i].points[j].Y());
            lines[j][i].SetY2(the_clusters[i].points[j].Y());
            n_lines[3]++;
          }
          if (the_clusters[i].plate_index == 8)
          {
            lines[j][i].SetY1(-1.92);
            lines[j][i].SetY2(1.92);
            lines[j][i].SetX1(the_clusters[i].points[j].X());
            lines[j][i].SetX2(the_clusters[i].points[j].X());
            n_lines[4]++;
          }
          if (the_clusters[i].plate_index == 9)
          {
            lines[j][i].SetX1(-3.84);
            lines[j][i].SetX2(3.84);
            lines[j][i].SetY1(the_clusters[i].points[j].Y());
            lines[j][i].SetY2(the_clusters[i].points[j].Y());
            n_lines[5]++;
          }
          if (the_clusters[i].plate_index == 10)
          {
            lines[j][i].SetY1(-3.84);
            lines[j][i].SetY2(3.84);
            lines[j][i].SetX1(the_clusters[i].points[j].X());
            lines[j][i].SetX2(the_clusters[i].points[j].X());
            n_lines[6]++;
          }
          if (the_clusters[i].plate_index == 11)
          {
            lines[j][i].SetX1(-3.84);
            lines[j][i].SetX2(3.84);
            lines[j][i].SetY1(the_clusters[i].points[j].Y());
            lines[j][i].SetY2(the_clusters[i].points[j].Y());
            n_lines[7]++;
          }
          if (the_clusters[i].plate_index == 12)
          {
            lines[j][i].SetY1(-3.84);
            lines[j][i].SetY2(3.84);
            lines[j][i].SetX1(the_clusters[i].points[j].X());
            lines[j][i].SetX2(the_clusters[i].points[j].X());
            n_lines[8]++;
          }
        }
      }
     
      for (int i = 0; i < n_track_options; i++)
      {
        h62->Fill(the_track_options[i].point_on_d_plate.X(), the_track_options[i].point_on_d_plate.Y());
      }
      
      break;
    }
  }
  
  fCanvas1->cd();
  h41->Draw("p");
  for (int i = 0; i < n_lines[0]; i++)
  {
    lines[i][0].Draw();
  }
  fCanvas1->Update();
  
  fCanvas2->cd();
  h51->Draw("p");
  for (int i = 0; i < n_lines[1]; i++)
  {
    lines[i][1].Draw();
  }
  fCanvas2->Update();
  
  fCanvas3->cd();
  h61->Draw("p");
  h62->Draw("same");
  for (int i = 0; i < n_lines[2]; i++)
  {
    lines[i][2].Draw();
  }
  fCanvas3->Update();
  
  fCanvas4->cd();
  h71->Draw("p");
  for (int i = 0; i < n_lines[3]; i++)
  {
    lines[i][3].Draw();
  }
  fCanvas4->Update();
  
  fCanvas5->cd();
  h81->Draw("p");
  for (int i = 0; i < n_lines[4]; i++)
  {
    lines[i][4].Draw();
  }
  fCanvas5->Update();
  
  fCanvas6->cd();
  h91->Draw("p");
  h93->Draw("same");
  for (int i = 0; i < n_lines[5]; i++)
  {
    lines[i][5].Draw();
  }
  fCanvas6->Update();
  
  fCanvas7->cd();
  h101->Draw("p");
  h103->Draw("same");
  for (int i = 0; i < n_lines[6]; i++)
  {
    lines[i][6].Draw();
  }
  fCanvas7->Update();
  
  fCanvas8->cd();
  h111->Draw("p");
  h113->Draw("same");
  for (int i = 0; i < n_lines[7]; i++)
  {
    lines[i][7].Draw();
  }
  fCanvas8->Update();
  
  fCanvas9->cd();
  h121->Draw("p");
  h123->Draw("same");
  for (int i = 0; i < n_lines[8]; i++)
  {
    lines[i][8].Draw();
  }
  fCanvas9->Update();
  
  if (h41) delete h41;
  if (h51) delete h51;
  if (h61) delete h61;
  if (h71) delete h71;
  if (h81) delete h81;
  if (h91) delete h91;
  if (h101) delete h101;
  if (h111) delete h111;
  if (h121) delete h121;
  if (h93) delete h93;
  if (h103) delete h103;
  if (h113) delete h113;
  if (h123) delete h123;
}
MyMainFrame::~MyMainFrame() {
  // Clean up used widgets: frames, buttons, layout hints
  fMain->Cleanup();
  in_file.close();
  delete fMain;
}
void example_2() {
  // Popup the GUI...
  new MyMainFrame(gClient->GetRoot(),2000,1600);
}
