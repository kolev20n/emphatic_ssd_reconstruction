void plot_vertex()
{
  double x1, y1, z1, x2, y2, z2;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 600, 400);
  
  TH1D* h1 = new TH1D("h1", "h1", 100, 0., 0.2);
  TH1D* h2 = new TH1D("h2", "h2", 100, 0., 0.2);
  TH1D* h3 = new TH1D("h3", "h3", 100, 0., 0.2);

  
  ifstream in_file("output_vertices.txt");
  
  while (in_file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2)
  {
    h1->Fill(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
  }
  
  in_file.close();
  
  in_file.open("output_vertices_misaligned.txt");
  
  while (in_file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2)
  {
    h2->Fill(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
  }
  
  in_file.close();

  in_file.open("output_vertices_highly_misaligned.txt");
  
  while (in_file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2)
  {
    h3->Fill(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
  }
  
  in_file.close();
  
  c1->cd();
  h1->Draw();
  h2->Draw("sames");
  h3->Draw("sames");

}
