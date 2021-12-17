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

void summary_plots()
{
  ifstream in_file;
  in_file.open("temp_debug_midstream_downstream_with_tracks_all.txt", ios::in);
  
  string card;
  int event_class;
  int cluster_multiplicity[13];
  double x, y, z;
  XYZVector mc_vertex, guessed_vertex, reco_vertex;
  int n_particles = 0;
  particle the_particles[100];
  int n_cluster_plates = 0;
  cluster the_clusters[100];
  int n_track_options;
  track_option the_track_options[1000];
  //int n_lines[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  //TLine lines[100][9];
  //const double sqrt2 = 0.707106781186548;
  //double t;
  int n_tracks;
  int n_track_hits[100];
  XYZVector the_track_hits[100][4];
  int n_track_hits_on_plate[4];
  
  int plate_index;
  
  double distance;
  
  int nev = 0;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 5, 26, 1200, 815);
  TH1F* h11 = new TH1F("h11", "h11", 100, 0, 0.05);
  TH1F* h12 = new TH1F("h12", "h12", 100, 0, 0.05);
  TH1F* h13 = new TH1F("h13", "h13", 100, 0, 0.05);
  TH1F* h14 = new TH1F("h14", "h14", 100, 0, 0.05);

  TCanvas *c2 = new TCanvas("c2", "c2", 5, 26, 1200, 815);
  TH1F* h21 = new TH1F("h21", "h21", 100, 0, 0.5);
  
  while (in_file >> card)
  {
    if (card == "event")
    {
      if (nev > 0)
      {
        //break;
      }
      else
      {
        nev++;
      }
      
      n_particles = 0;
      n_cluster_plates = 0;
      n_track_options = 0;
      n_tracks = 0;
      n_track_hits_on_plate[0] = 0;
      n_track_hits_on_plate[1] = 0;
      n_track_hits_on_plate[2] = 0;
      n_track_hits_on_plate[3] = 0;
      
      for (int i = 0; i < 100; i++)
      {
        the_particles[i].n_plates = 0;
        the_clusters[i].n_clusters = 0;
        n_track_hits[i] = 0;
      }
      
      in_file >> event_class;
      for (int i = 0; i < 13; i++)
      {
        in_file >> cluster_multiplicity[i];
      }
      in_file >> x >> y >> z;
      mc_vertex.SetCoordinates(x, y, z);
      in_file >> x >> y >> z;
      guessed_vertex.SetCoordinates(x, y, z);
      in_file >> x >> y >> z;
      reco_vertex.SetCoordinates(x, y, z);
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
        in_file >> n_track_hits[i];
        for (int j = 0; j < n_track_hits[i]; j++)
        {
          in_file >> plate_index >> x >> y >> z;
          if (plate_index == 9)
          {
            the_track_hits[n_track_hits_on_plate[0]][0].SetCoordinates(x, y, z);
            n_track_hits_on_plate[0]++;
          }
          else if (plate_index == 10)
          {
            the_track_hits[n_track_hits_on_plate[1]][1].SetCoordinates(x, y, z);
            n_track_hits_on_plate[1]++;
          }
          else if (plate_index == 11)
          {
            the_track_hits[n_track_hits_on_plate[2]][2].SetCoordinates(x, y, z);
            n_track_hits_on_plate[2]++;
          }
          else if (plate_index == 12)
          {
            the_track_hits[n_track_hits_on_plate[3]][3].SetCoordinates(x, y, z);
            n_track_hits_on_plate[3]++;
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
          if (the_particles[i].plate_index[j] == 9)
          {
            for (int k = 0; k < n_track_hits_on_plate[0]; k++)
            {
              distance = sqrt((the_track_hits[k][0].X() - the_particles[i].plate_hits[j].X()) *(the_track_hits[k][0].X() - the_particles[i].plate_hits[j].X()) +
                              (the_track_hits[k][0].Y() - the_particles[i].plate_hits[j].Y()) *(the_track_hits[k][0].Y() - the_particles[i].plate_hits[j].Y()));
              if (distance < 0.5) h11->Fill(distance);
            }
          }
          if (the_particles[i].plate_index[j] == 10)
          {
            for (int k = 0; k < n_track_hits_on_plate[1]; k++)
            {
              distance = sqrt((the_track_hits[k][1].X() - the_particles[i].plate_hits[j].X()) * (the_track_hits[k][1].X() - the_particles[i].plate_hits[j].X()) +
                              (the_track_hits[k][1].Y() - the_particles[i].plate_hits[j].Y()) * (the_track_hits[k][1].Y() - the_particles[i].plate_hits[j].Y()));
              if (distance < 0.5) h12->Fill(distance);
            }
          }
          if (the_particles[i].plate_index[j] == 11)
          {
            for (int k = 0; k < n_track_hits_on_plate[2]; k++)
            {
              distance = sqrt((the_track_hits[k][2].X() - the_particles[i].plate_hits[j].X()) * (the_track_hits[k][2].X() - the_particles[i].plate_hits[j].X()) +
                              (the_track_hits[k][2].Y() - the_particles[i].plate_hits[j].Y()) * (the_track_hits[k][2].Y() - the_particles[i].plate_hits[j].Y()));
              if (distance < 0.5) h13->Fill(distance);
            }
          }
          if (the_particles[i].plate_index[j] == 12)
          {
            for (int k = 0; k < n_track_hits_on_plate[3]; k++)
            {
              distance = sqrt((the_track_hits[k][3].X() - the_particles[i].plate_hits[j].X()) * (the_track_hits[k][3].X() - the_particles[i].plate_hits[j].X()) +
                              (the_track_hits[k][3].Y() - the_particles[i].plate_hits[j].Y()) * (the_track_hits[k][3].Y() - the_particles[i].plate_hits[j].Y()));
              if (distance < 0.5) h14->Fill(distance);
            }
          }
        }
      }
      
      if (reco_vertex.X() == 0. && reco_vertex.Y() == 0. && reco_vertex.Z() == 0.)
      {
      }
      else
      {
cout << sqrt((reco_vertex - mc_vertex).Mag2()) << endl;
        h21->Fill(sqrt((reco_vertex - mc_vertex).Mag2()));
      }
    }
  }
  
  c1->cd();
  h11->Draw("");
  h12->Draw("sames");
  h13->Draw("sames");
  h14->Draw("sames");

  c2->cd();
  h21->Draw("");
  
  return;
}
