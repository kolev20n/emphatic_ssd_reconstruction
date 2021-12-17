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
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

using namespace std;
using namespace ROOT::Math;

class particle
{
public:
  int id;
  int code;
  double ekin_vertex;
  int n_plates;
  int plate_index[16];
  XYZVector plate_hits[16];
};

class cluster
{
public:
  int plate_index;
  int n_clusters;
  XYZVector points[100];
  XYZVector directions[100];
};

class MyMainFrame
{
  RQ_OBJECT("MyMainFrame")
private:
  TGMainFrame         *fMain, *my_frame;
  // TRootEmbeddedCanvas *ec_0, *ec_1, *ec_2, *ec_3, *ec_4, *ec_5, *ec_6, *ec_7, *ec_8, *ec_9, *ec_10, *ec_11;
  TRootEmbeddedCanvas* ec[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  ifstream in_file;
  
public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();
  void DoDraw();
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h)
{
  string canvas_names[12] = {"mid_y1", "mid_x1", "mid_d2", "mid_d1", "mid_y2", "mid_x2", "down_y1", "down_x1", "down_d2", "down_d1", "down_y2", "down_x2"};
  
  in_file.open("output_debug_22_4.txt", ios::in);
  
  fMain = new TGMainFrame(p,w,h);
  my_frame = new TGMainFrame(p,w,h);
  
  TGVerticalFrame *f3 = new TGVerticalFrame(fMain, 1200, 2000);
  TGHorizontalFrame *f1 = new TGHorizontalFrame(f3, 1200, 800);
  TGHorizontalFrame *f2 = new TGHorizontalFrame(f3, 1200, 800);

  /*
  ec_0 = new TRootEmbeddedCanvas(canvas_names[0].c_str(), f1, 1200, 600);
  ec_1 = new TRootEmbeddedCanvas(canvas_names[1].c_str(), f1, 1200, 600);
  ec_2 = new TRootEmbeddedCanvas(canvas_names[2].c_str(), f1, 1200, 600);
  ec_3 = new TRootEmbeddedCanvas(canvas_names[3].c_str(), f1, 1200, 600);
  ec_4 = new TRootEmbeddedCanvas(canvas_names[4].c_str(), f1, 1200, 600);
  ec_5 = new TRootEmbeddedCanvas(canvas_names[5].c_str(), f1, 1200, 600);

  ec_6 = new TRootEmbeddedCanvas(canvas_names[6].c_str(), f2, 1200, 600);
  ec_7 = new TRootEmbeddedCanvas(canvas_names[7].c_str(), f2, 1200, 600);
  ec_8 = new TRootEmbeddedCanvas(canvas_names[8].c_str(), f2, 1200, 600);
  ec_9 = new TRootEmbeddedCanvas(canvas_names[9].c_str(), f2, 1200, 600);
  ec_10 = new TRootEmbeddedCanvas(canvas_names[10].c_str(), f2, 1200, 600);
  ec_11 = new TRootEmbeddedCanvas(canvas_names[11].c_str(), f2, 1200, 600);
  */

  
  for (int i = 0; i < 12; i++)
  {
    if (i < 6)
    {
      ec[i] = new TRootEmbeddedCanvas(canvas_names[i].c_str(), f1, 1200, 600);
      f1->AddFrame(ec[i], new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));
    }
    else
    {
      ec[i] = new TRootEmbeddedCanvas(canvas_names[i].c_str(), f2, 1200, 600);
      f2->AddFrame(ec[i], new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));
    }
  }
  
  f3->AddFrame(f1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  f3->AddFrame(f2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  fMain->AddFrame(f3, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain, 200, 40);
  TGTextButton *draw = new TGTextButton(hframe, "&Draw");
  draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  TGTextButton *exit = new TGTextButton(hframe, "&Exit", "gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
  
  fMain->SetWindowName("Simple Example");
  my_frame->SetWindowName("Frame");
  
  fMain->MapSubwindows();
  
  fMain->Resize(fMain->GetDefaultSize());
  
  fMain->MapWindow();
  my_frame->MapWindow();
}

void MyMainFrame::DoDraw()
{
  string card;
  int event_class;
  int cluster_multiplicity[16];
  double x, y, z;
  XYZVector mc_vertex, guessed_vertex;
  int n_particles = 0;
  particle the_particles[100];
  int n_cluster_plates = 0;
  cluster the_clusters[100];
  int n_track_options;
  int n_lines[12] = {0, 0, 0, 0, 0};
  TLine lines[100][12];
  const double sqrt2 = 0.707106781186548;
  const double to_radians = 3.141592 / 180.;
  double t;
  int dummy_int;
  
  double pattern_quality;
  int n_tracks;
  XYZVector tracks[100];
  
  RotationZ rz1, rz2;
  rz1.SetAngle( 45. * to_radians);
  rz2.SetAngle(-45. * to_radians);
  XYZVector v1, v2;
  
  string histogram_names[12] = {"plate_mid_y1", "plate_mid_x1", "plate_mid_d2", "plate_mid_d1", "plate_mid_y2", "plate_mid_x2", "plate_down_y1", "plate_down_x1", "plate_down_d2", "plate_down_d1", "plate_down_y2", "plate_down_x2"};
  string track_histogram_names[12] = {"track_mid_y1", "track_mid_x1", "track_mid_d2", "track_mid_d1", "track_mid_y2", "track_mid_x2", "track_down_y1", "track_down_x1", "track_down_d2", "track_down_d1", "track_down_y2", "track_down_x2"};

  TH2D* hp[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  TH2D* ht[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  for (int i = 0; i < 12; i++)
  {
    if (i < 6)
    {
      hp[i] = new TH2D(histogram_names[i].c_str(), histogram_names[i].c_str(),   640, -1.92, 1.92, 640,  -1.92, 1.92);
      ht[i] = new TH2D(track_histogram_names[i].c_str(), track_histogram_names[i].c_str(),   640, -1.92, 1.92, 640,  -1.92, 1.92);
    }
    else
    {
      hp[i] = new TH2D(histogram_names[i].c_str(), histogram_names[i].c_str(),   640, -3.84, 3.84, 640,  -3.84, 3.84);
      ht[i] = new TH2D(track_histogram_names[i].c_str(), track_histogram_names[i].c_str(),   640, -3.84, 3.84, 640,  -3.84, 3.84);
    }
  }
  
  for (int i = 0; i < 12; i++)
  {
    hp[i]->SetMarkerColor(4);
    hp[i]->SetMarkerStyle(20);
    hp[i]->SetMarkerSize(1);
    hp[i]->SetStats(FALSE);
  }
  
  for (int i = 0; i < 12; i++)
  {
    ht[i]->SetMarkerColor(2);
    ht[i]->SetMarkerStyle(43);
    ht[i]->SetMarkerSize(2);
    ht[i]->SetStats(FALSE);
  }

  TCanvas* the_canvases[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  for (int i = 0; i < 12; i++)
  {
    the_canvases[i] = ec[i]->GetCanvas();
  }

  while (in_file >> card)
  {
//cout << card << endl;
    if (card == "event")
    {
      in_file >> card >> event_class;
    }
    else if (card == "particle")
    {
      in_file >> the_particles[n_particles].id;
      in_file >> the_particles[n_particles].code;
      in_file >> the_particles[n_particles].ekin_vertex;
      in_file >> the_particles[n_particles].n_plates;
      for (int i = 0; i < the_particles[n_particles].n_plates; i++)
      {
        in_file >> dummy_int;

        if (dummy_int < 4)
        {
          in_file >> x >> y >> z;
        }
        else if (dummy_int < 10)
        {
          the_particles[n_particles].plate_index[i] = dummy_int;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 10 || dummy_int == 11)
        {
          the_particles[n_particles].plate_index[i] = 10;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 12 || dummy_int == 13)
        {
          the_particles[n_particles].plate_index[i] = 11;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 14 || dummy_int == 15)
        {
          the_particles[n_particles].plate_index[i] = 12;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 16 || dummy_int == 17)
        {
          the_particles[n_particles].plate_index[i] = 13;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 18 || dummy_int == 19)
        {
          the_particles[n_particles].plate_index[i] = 14;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
        else if (dummy_int == 20 || dummy_int == 21)
        {
          the_particles[n_particles].plate_index[i] = 15;
          in_file >> x >> y >> z;
          the_particles[n_particles].plate_hits[i].SetCoordinates(x, y, z);
        }
      }
      if (dummy_int > 3)
      {
        n_particles++;
      }
    }
    else if (card == "number_of_clusters")
    {
      for (int i = 0; i < 16; i++)
      {
        in_file >> cluster_multiplicity[i];
      }
    }
    else if (card == "mc_vertex")
    {
      in_file >> x >> y >> z;
      mc_vertex.SetCoordinates(x, y, z);
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
    else if (card == "pattern")
    {
      in_file >> n_tracks >> pattern_quality;
    }
    else if (card == "vertex")
    {
      in_file >> x >> y >> z;
    }
    else if (card == "tracks")
    {
      for (int i = 0; i < n_tracks; i++)
      {
        in_file >> dummy_int >> x >> y >> z;
        tracks[i].SetCoordinates(x, y, z);
      }
      for (int i = 0; i < n_tracks; i++)
      {
        for (int j = 0; j < 16; j++)
        {
          in_file >> x >> y >> z;
//cout << x << y << z << endl;
//cout << "+++" << endl;
          if (!(x == 0. && y == 0. && z == 0.) && j > 3)
          {
//cout << x << y << z << endl;
            ht[j - 4]->Fill(x, y);
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
          hp[the_particles[i].plate_index[j] - 4]->Fill(the_particles[i].plate_hits[j].X(), the_particles[i].plate_hits[j].Y());
        }
      }
      
      for (int i = 4; i < n_cluster_plates; i++)
      {
        for (int j = 0; j < the_clusters[i].n_clusters; j++)
        {
          if (the_clusters[i].plate_index == 4)
          {
            lines[j][i - 4].SetX1(-1.92);
            lines[j][i - 4].SetX2(1.92);
            lines[j][i - 4].SetY1(the_clusters[i].points[j].Y());
            lines[j][i - 4].SetY2(the_clusters[i].points[j].Y());
            n_lines[0]++;
          }
          else if (the_clusters[i].plate_index == 5)
          {
            lines[j][i - 4].SetY1(-1.92);
            lines[j][i - 4].SetY2(1.92);
            lines[j][i - 4].SetX1(the_clusters[i].points[j].X());
            lines[j][i - 4].SetX2(the_clusters[i].points[j].X());
            n_lines[1]++;
          }
          else if (the_clusters[i].plate_index == 6)
          {
            if (the_clusters[i].points[j].X() < 0.)
            {
              t = (-1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetX1(-1.92);
              lines[j][i - 4].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
            
              t = (-1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetY2(-1.92);
              n_lines[2]++;
            }
            else
            {
              t = (1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetX1(1.92);
              lines[j][i - 4].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetY2(1.92);
              n_lines[2]++;
            }
          }
          else if (the_clusters[i].plate_index == 7)
          {
            if (the_clusters[i].points[j].X() > 0.)
            {
              t = (1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetY1(1.92);
              lines[j][i - 4].SetX1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (-1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetY2((the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t));
              lines[j][i - 4].SetX2(-1.92);
              n_lines[3]++;
            }
            else
            {
              t = (-1.92 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetY1(-1.92);
              lines[j][i - 4].SetX1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (1.92 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetY2((the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t));
              lines[j][i - 4].SetX2(1.92);
              n_lines[3]++;
            }
          }
          else if (the_clusters[i].plate_index == 8)
          {
            lines[j][i - 4].SetX1(-1.92);
            lines[j][i - 4].SetX2(1.92);
            lines[j][i - 4].SetY1(the_clusters[i].points[j].Y());
            lines[j][i - 4].SetY2(the_clusters[i].points[j].Y());
            n_lines[4]++;
          }
          else if (the_clusters[i].plate_index == 9)
          {
            lines[j][i - 4].SetY1(-1.92);
            lines[j][i - 4].SetY2(1.92);
            lines[j][i - 4].SetX1(the_clusters[i].points[j].X());
            lines[j][i - 4].SetX2(the_clusters[i].points[j].X());
            n_lines[5]++;
          }
          else if (the_clusters[i].plate_index == 10)
          {
            lines[j][i - 4].SetX1(-3.84);
            lines[j][i - 4].SetX2(3.84);
            lines[j][i - 4].SetY1(the_clusters[i].points[j].Y());
            lines[j][i - 4].SetY2(the_clusters[i].points[j].Y());
            n_lines[6]++;
          }
          else if (the_clusters[i].plate_index == 11)
          {
            lines[j][i - 4].SetY1(-3.84);
            lines[j][i - 4].SetY2(3.84);
            lines[j][i - 4].SetX1(the_clusters[i].points[j].X());
            lines[j][i - 4].SetX2(the_clusters[i].points[j].X());
            n_lines[7]++;
          }
          else if (the_clusters[i].plate_index == 12)
          {
            if (the_clusters[i].points[j].X() < 0.)
            {
              t = (-3.84 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetX1(-3.84);
              lines[j][i - 4].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);

              t = (-3.84 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetY2(-3.84);
              n_lines[8]++;
            }
            else
            {
              t = (3.84 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetX1(3.84);
              lines[j][i - 4].SetY1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (3.84 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetX2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetY2(3.84);
              n_lines[8]++;
            }
          }
          else if (the_clusters[i].plate_index == 13)
          {
            if (the_clusters[i].points[j].X() > 0.)
            {
              t = (3.84 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetY1(3.84);
              lines[j][i - 4].SetX1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (-3.84 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetY2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetX2(-3.84);
              n_lines[9]++;
            }
            else
            {
              t = (-3.84 - the_clusters[i].points[j].X()) / the_clusters[i].directions[j].X();
              
              lines[j][i - 4].SetY1(-3.84);
              lines[j][i - 4].SetX1(the_clusters[i].points[j].Y() + the_clusters[i].directions[j].Y() * t);
              
              t = (3.84 - the_clusters[i].points[j].Y()) / the_clusters[i].directions[j].Y();
              lines[j][i - 4].SetY2(the_clusters[i].points[j].X() + the_clusters[i].directions[j].X() * t);
              lines[j][i - 4].SetX2(3.84);
              n_lines[9]++;
            }
          }
          else if (the_clusters[i].plate_index == 14)
          {
            lines[j][i - 4].SetX1(-3.84);
            lines[j][i - 4].SetX2(3.84);
            lines[j][i - 4].SetY1(the_clusters[i].points[j].Y());
            lines[j][i - 4].SetY2(the_clusters[i].points[j].Y());
            n_lines[10]++;
          }
          else if (the_clusters[i].plate_index == 15)
          {
            lines[j][i - 4].SetY1(-3.84);
            lines[j][i - 4].SetY2(3.84);
            lines[j][i - 4].SetX1(the_clusters[i].points[j].X());
            lines[j][i - 4].SetX2(the_clusters[i].points[j].X());
            n_lines[11]++;
          }
        }
      }

      
      break;
    }
  }
  
  for (int i = 0; i < 12; i++)
  {
    the_canvases[i]->cd();
    hp[i]->Draw("p");
    for (int j = 0; j < n_lines[i]; j++)
    {
      lines[j][i].Draw();
    }
    ht[i]->Draw("same");
    the_canvases[i]->Update();
  }
  
  for (int i = 0; i < 12; i++)
  {
    if (hp[i]) delete hp[i];
    if (ht[i]) delete ht[i];
  }

}

MyMainFrame::~MyMainFrame()
{
  fMain->Cleanup();
  in_file.close();
  delete fMain;
}

void ssd_display_22_4()
{
  new MyMainFrame(gClient->GetRoot(), 2000, 1600);
}
