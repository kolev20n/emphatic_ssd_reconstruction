#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>
#include <TColor.h>
#include <TRandom.h>
#include <TMath.h>
#include <TEveBox.h>
#include <TEveArrow.h>
#include <TEveBrowser.h>

#include <string>
#include <vector>

using namespace std;

class input_root_data_structure
{
public:
  input_root_data_structure();
  ~input_root_data_structure();
  
  std::vector<int> *plate_number;
  std::vector<int> *strip_number;
  std::vector<double> *total_energy_values;
  std::vector<double> *non_ionization_energy_values;
  std::vector<std::vector<double>> *contributing_tracks;
  
  std::vector<int> *track_id;
  std::vector<int> *particle_code;
  std::vector<int> *parent_track_id;
  std::vector<std::string> *creator_process;
  std::vector<double> *x_vertex;
  std::vector<double> *y_vertex;
  std::vector<double> *z_vertex;
  std::vector<double> *px_vertex;
  std::vector<double> *py_vertex;
  std::vector<double> *pz_vertex;
  std::vector<double> *ekin_vertex;
  std::vector<std::vector<double>> *visited_plates;
  
  std::vector<int> *step_track_id;
  std::vector<std::string> *step_volume_name;
  std::vector<int> *step_plate_id;
  std::vector<double> *step_x_i;
  std::vector<double> *step_y_i;
  std::vector<double> *step_z_i;
  std::vector<double> *step_px_i;
  std::vector<double> *step_py_i;
  std::vector<double> *step_pz_i;
  std::vector<double> *step_x_f;
  std::vector<double> *step_y_f;
  std::vector<double> *step_z_f;
  std::vector<double> *step_px_f;
  std::vector<double> *step_py_f;
  std::vector<double> *step_pz_f;
  std::vector<double> *step_initial_energy;
  std::vector<double> *step_final_energy;
  std::vector<std::string> *step_process_name;
};

input_root_data_structure::input_root_data_structure()
:plate_number(0), strip_number(0), total_energy_values(0), non_ionization_energy_values(0), contributing_tracks(0),
track_id(0), particle_code(0),
parent_track_id(0), creator_process(0), x_vertex(0),
y_vertex(0), z_vertex(0), px_vertex(0),
py_vertex(0), pz_vertex(0), ekin_vertex(0),
visited_plates(0), step_track_id(0), step_volume_name(0),
step_plate_id(0), step_x_i(0), step_y_i(0), step_z_i(0),
step_px_i(0), step_py_i(0), step_pz_i(0),
step_x_f(0), step_y_f(0), step_z_f(0),
step_px_f(0), step_py_f(0), step_pz_f(0),
step_initial_energy(0), step_final_energy(0), step_process_name(0)
{
  cout << "Created the input root data structure." << endl;
}

int my_plot(input_root_data_structure* the_input_root_data_structure, int current_event)
{
  /*
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
  */
  
  double um = 1.e-4;
  
  double midstream_half_size = 1.92;
  double downstream_half_size = 3.84;

  double midstream_shift = 6.0;
  double downstream_shift = 10.0;

  int midstream_number_of_strips = 640;
  int downstream_number_of_strips = 1280;
  
  double dummy_double;
  
  TEveBox* m0 = new TEveBox();
  m0->SetVertex(0, -midstream_half_size, -midstream_half_size, 0.);
  m0->SetVertex(1, -midstream_half_size, midstream_half_size, 0.);
  m0->SetVertex(2, -midstream_half_size, midstream_half_size, 0.02);
  m0->SetVertex(3, -midstream_half_size, -midstream_half_size, 0.02);
  m0->SetVertex(4, midstream_half_size, -midstream_half_size, 0.);
  m0->SetVertex(5, midstream_half_size, midstream_half_size, 0.);
  m0->SetVertex(6, midstream_half_size, midstream_half_size, 0.02);
  m0->SetVertex(7, midstream_half_size, -midstream_half_size, 0.02);
  
  TEveBox* m1 = new TEveBox();
  m1->SetVertex(0, -midstream_half_size + midstream_shift, -midstream_half_size, 0.);
  m1->SetVertex(1, -midstream_half_size + midstream_shift, midstream_half_size, 0.);
  m1->SetVertex(2, -midstream_half_size + midstream_shift, midstream_half_size, 0.02);
  m1->SetVertex(3, -midstream_half_size + midstream_shift, -midstream_half_size, 0.02);
  m1->SetVertex(4, midstream_half_size + midstream_shift, -midstream_half_size, 0.);
  m1->SetVertex(5, midstream_half_size + midstream_shift, midstream_half_size, 0.);
  m1->SetVertex(6, midstream_half_size + midstream_shift, midstream_half_size, 0.02);
  m1->SetVertex(7, midstream_half_size + midstream_shift, -midstream_half_size, 0.02);
  
  TEveBox* m2 = new TEveBox();
  m2->SetVertex(0, -midstream_half_size + 2.*midstream_shift, -midstream_half_size, 0.);
  m2->SetVertex(1, -midstream_half_size + 2.*midstream_shift, midstream_half_size, 0.);
  m2->SetVertex(2, -midstream_half_size + 2.*midstream_shift, midstream_half_size, 0.02);
  m2->SetVertex(3, -midstream_half_size + 2.*midstream_shift, -midstream_half_size, 0.02);
  m2->SetVertex(4, midstream_half_size + 2.*midstream_shift, -midstream_half_size, 0.);
  m2->SetVertex(5, midstream_half_size + 2.*midstream_shift, midstream_half_size, 0.);
  m2->SetVertex(6, midstream_half_size + 2.*midstream_shift, midstream_half_size, 0.02);
  m2->SetVertex(7, midstream_half_size + 2.*midstream_shift, -midstream_half_size, 0.02);
  
  TEveBox* m3 = new TEveBox();
  m3->SetVertex(0, -midstream_half_size + 3.*midstream_shift, -midstream_half_size, 0.);
  m3->SetVertex(1, -midstream_half_size + 3.*midstream_shift, midstream_half_size, 0.);
  m3->SetVertex(2, -midstream_half_size + 3.*midstream_shift, midstream_half_size, 0.02);
  m3->SetVertex(3, -midstream_half_size + 3.*midstream_shift, -midstream_half_size, 0.02);
  m3->SetVertex(4, midstream_half_size + 3.*midstream_shift, -midstream_half_size, 0.);
  m3->SetVertex(5, midstream_half_size + 3.*midstream_shift, midstream_half_size, 0.);
  m3->SetVertex(6, midstream_half_size + 3.*midstream_shift, midstream_half_size, 0.02);
  m3->SetVertex(7, midstream_half_size + 3.*midstream_shift, -midstream_half_size, 0.02);
  
  TEveBox* m4 = new TEveBox();
  m4->SetVertex(0, -midstream_half_size + 4.*midstream_shift, -midstream_half_size, 0.);
  m4->SetVertex(1, -midstream_half_size + 4.*midstream_shift, midstream_half_size, 0.);
  m4->SetVertex(2, -midstream_half_size + 4.*midstream_shift, midstream_half_size, 0.02);
  m4->SetVertex(3, -midstream_half_size + 4.*midstream_shift, -midstream_half_size, 0.02);
  m4->SetVertex(4, midstream_half_size + 4.*midstream_shift, -midstream_half_size, 0.);
  m4->SetVertex(5, midstream_half_size + 4.*midstream_shift, midstream_half_size, 0.);
  m4->SetVertex(6, midstream_half_size + 4.*midstream_shift, midstream_half_size, 0.02);
  m4->SetVertex(7, midstream_half_size + 4.*midstream_shift, -midstream_half_size, 0.02);
  
  TEveBox* m5 = new TEveBox();
  m5->SetVertex(0, -midstream_half_size + 5.*midstream_shift, -midstream_half_size, 0.);
  m5->SetVertex(1, -midstream_half_size + 5.*midstream_shift, midstream_half_size, 0.);
  m5->SetVertex(2, -midstream_half_size + 5.*midstream_shift, midstream_half_size, 0.02);
  m5->SetVertex(3, -midstream_half_size + 5.*midstream_shift, -midstream_half_size, 0.02);
  m5->SetVertex(4, midstream_half_size + 5.*midstream_shift, -midstream_half_size, 0.);
  m5->SetVertex(5, midstream_half_size + 5.*midstream_shift, midstream_half_size, 0.);
  m5->SetVertex(6, midstream_half_size + 5.*midstream_shift, midstream_half_size, 0.02);
  m5->SetVertex(7, midstream_half_size + 5.*midstream_shift, -midstream_half_size, 0.02);
  
  TEvePointSet* ps = new TEvePointSet();
  
  //cout << the_input_root_data_structure->plate_number->size() << endl;
  
  for (int i = 0; i < the_input_root_data_structure->plate_number->size(); i++)
  {
    if (the_input_root_data_structure->plate_number->at(i) == 4)
    {
      //cout << "!!" << endl;
      dummy_double = (the_input_root_data_structure->strip_number->at(i) - (((double)(midstream_number_of_strips)) - 1.) / 2.) * 60 * um;
      ps->SetNextPoint(0., dummy_double, 0.);
    }
    else if (the_input_root_data_structure->plate_number->at(i) == 6)
    {
      dummy_double = (the_input_root_data_structure->strip_number->at(i) - (((double)(midstream_number_of_strips)) - 1.) / 2.) * 60 * um;
      ps->SetNextPoint(0. + 2.*midstream_shift, dummy_double, 0.);
    }
    else if (the_input_root_data_structure->plate_number->at(i) == 8)
    {
      dummy_double = (the_input_root_data_structure->strip_number->at(i) - (((double)(midstream_number_of_strips)) - 1.) / 2.) * 60 * um;
      ps->SetNextPoint(0. + 4.*midstream_shift, dummy_double, 0.);
    }
  }
  
  ps->SetMarkerColor(3);
  ps->SetMarkerSize(2);
  ps->SetMarkerStyle(4);
  
  
  if (gEve->GetViewers()) gEve->GetViewers()->DeleteAnnotations();
  if (gEve->GetCurrentEvent()) gEve->GetCurrentEvent()->DestroyElements();
  
  gEve->AddElement(m0);
  gEve->AddElement(m1);
  gEve->AddElement(m2);
  gEve->AddElement(m3);
  gEve->AddElement(m4);
  gEve->AddElement(m5);

  gEve->AddElement(ps);
  
  gEve->Redraw3D(true);
  
  //app->Terminate(0);
  
  return 0;
}

class EvNavHandler
{
public:
  EvNavHandler(input_root_data_structure* an_input_root_data_structure, TTree* an_input_ssd_energy_tree, TTree* an_input_particle_truth_tree, TTree* an_input_steps_truth_tree)
  {
    the_input_root_data_structure = an_input_root_data_structure;
    the_input_ssd_energy_tree = an_input_ssd_energy_tree;
    the_input_particle_truth_tree = an_input_particle_truth_tree;
    the_input_steps_truth_tree = an_input_steps_truth_tree;
    current_event = -1;
    total_events = the_input_ssd_energy_tree->GetEntries();
  }
  void Fwd()
  {
    current_event++;
    if (current_event < total_events)
    {
      the_input_ssd_energy_tree->GetEntry(current_event);
      the_input_particle_truth_tree->GetEntry(current_event);
      the_input_steps_truth_tree->GetEntry(current_event);
    }
    else
    {
      current_event--;
      cout << "Already at last event..." << endl;
    }
    
    //cout << current_event << endl;
    //cout << the_input_root_data_structure->plate_number->size() << endl;
    my_plot(the_input_root_data_structure, current_event);
  }
  void Bck()
  {
    current_event--;
    if (current_event >= 0)
    {
      the_input_ssd_energy_tree->GetEntry(current_event);
      the_input_particle_truth_tree->GetEntry(current_event);
      the_input_steps_truth_tree->GetEntry(current_event);
    }
    else
    {
      current_event++;
      cout << "Already at first event..." << endl;
    }
    
    //cout << current_event << endl;
    //cout << the_input_root_data_structure->plate_number->size() << endl;
    my_plot(the_input_root_data_structure, current_event);
  }
private:
  input_root_data_structure* the_input_root_data_structure;
  TTree* the_input_ssd_energy_tree;
  TTree* the_input_particle_truth_tree;
  TTree* the_input_steps_truth_tree;
  int current_event;
  int total_events;
};

void make_gui(input_root_data_structure* the_input_root_data_structure, TTree* the_input_ssd_energy_tree, TTree* the_input_particle_truth_tree, TTree* the_input_steps_truth_tree)
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
    EvNavHandler    *fh = new EvNavHandler(the_input_root_data_structure, the_input_ssd_energy_tree, the_input_particle_truth_tree, the_input_steps_truth_tree);
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


void run_vis()
{
  int current_event = 0;
  int total_events;
  
  input_root_data_structure* the_input_root_data_structure = new input_root_data_structure();
  
  TFile* in_file = new TFile("../../data/production_30gev_proton_carbon_real_1_single_track_small.root", "READ");
  TTree* the_input_ssd_energy_tree = (TTree*) in_file->Get("ssd_energy");
  TTree* the_input_particle_truth_tree = (TTree*) in_file->Get("particle_truth");
  TTree* the_input_steps_truth_tree = (TTree*) in_file->Get("steps_truth");
  
  the_input_ssd_energy_tree->SetBranchAddress("plate_id", &the_input_root_data_structure->plate_number);
  the_input_ssd_energy_tree->SetBranchAddress("strip_id", &the_input_root_data_structure->strip_number);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_total_energy_deposited", &the_input_root_data_structure->total_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_non_ionization_energy_deposited", &the_input_root_data_structure->non_ionization_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("contributing_tracks", &the_input_root_data_structure->contributing_tracks);
  
  the_input_particle_truth_tree->SetBranchAddress("track_id", &the_input_root_data_structure->track_id);
  the_input_particle_truth_tree->SetBranchAddress("particle_code", &the_input_root_data_structure->particle_code);
  the_input_particle_truth_tree->SetBranchAddress("parent_track_id", &the_input_root_data_structure->parent_track_id);
  the_input_particle_truth_tree->SetBranchAddress("creator_process", &the_input_root_data_structure->creator_process);
  the_input_particle_truth_tree->SetBranchAddress("x_vertex", &the_input_root_data_structure->x_vertex);
  the_input_particle_truth_tree->SetBranchAddress("y_vertex", &the_input_root_data_structure->y_vertex);
  the_input_particle_truth_tree->SetBranchAddress("z_vertex", &the_input_root_data_structure->z_vertex);
  the_input_particle_truth_tree->SetBranchAddress("px_vertex", &the_input_root_data_structure->px_vertex);
  the_input_particle_truth_tree->SetBranchAddress("py_vertex", &the_input_root_data_structure->py_vertex);
  the_input_particle_truth_tree->SetBranchAddress("pz_vertex", &the_input_root_data_structure->pz_vertex);
  the_input_particle_truth_tree->SetBranchAddress("ekin_vertex", &the_input_root_data_structure->ekin_vertex);
  the_input_particle_truth_tree->SetBranchAddress("visited_plates", &the_input_root_data_structure->visited_plates);
  
  the_input_steps_truth_tree->SetBranchAddress("step_track_id", &the_input_root_data_structure->step_track_id);
  the_input_steps_truth_tree->SetBranchAddress("step_volume_name", &the_input_root_data_structure->step_volume_name);
  the_input_steps_truth_tree->SetBranchAddress("step_plate_id", &the_input_root_data_structure->step_plate_id);
  the_input_steps_truth_tree->SetBranchAddress("step_x_i", &the_input_root_data_structure->step_x_i);
  the_input_steps_truth_tree->SetBranchAddress("step_y_i", &the_input_root_data_structure->step_y_i);
  the_input_steps_truth_tree->SetBranchAddress("step_z_i", &the_input_root_data_structure->step_z_i);
  the_input_steps_truth_tree->SetBranchAddress("step_px_i", &the_input_root_data_structure->step_px_i);
  the_input_steps_truth_tree->SetBranchAddress("step_py_i", &the_input_root_data_structure->step_py_i);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_i", &the_input_root_data_structure->step_pz_i);
  the_input_steps_truth_tree->SetBranchAddress("step_x_f", &the_input_root_data_structure->step_x_f);
  the_input_steps_truth_tree->SetBranchAddress("step_y_f", &the_input_root_data_structure->step_y_f);
  the_input_steps_truth_tree->SetBranchAddress("step_z_f", &the_input_root_data_structure->step_z_f);
  the_input_steps_truth_tree->SetBranchAddress("step_px_f", &the_input_root_data_structure->step_px_f);
  the_input_steps_truth_tree->SetBranchAddress("step_py_f", &the_input_root_data_structure->step_py_f);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_f", &the_input_root_data_structure->step_pz_f);
  the_input_steps_truth_tree->SetBranchAddress("step_initial_energy", &the_input_root_data_structure->step_initial_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_final_energy", &the_input_root_data_structure->step_final_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_process_name", &the_input_root_data_structure->step_process_name);
  
  total_events = the_input_ssd_energy_tree->GetEntries();
  cout << "Total events: " << total_events << endl;
  
  TEveManager::Create();

  gEve->GetBrowser()->GetTabRight()->SetTab(1);
  make_gui(the_input_root_data_structure, the_input_ssd_energy_tree, the_input_particle_truth_tree, the_input_steps_truth_tree);
  
  

}
