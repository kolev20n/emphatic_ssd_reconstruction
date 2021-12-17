#include <iostream>
#include <utility>
#include <algorithm>
#include <math.h>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include "track_finding_manager.h"
#include "partial_geometry.h"
#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"
#include "my_straight_line_fcn.h"
#include "my_straight_line_fcn_with_point.h"
#include "my_global_fit_no_field_fcn.h"

using namespace std;
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

track_finding_manager::track_finding_manager()
:the_geometry(0), the_reconstruction_algorithm(0)
{
  cout << "Created the_track_finding_manager." << endl;

  the_geometry = the_global_data_dispatcher->get_partial_geometry();
  
  the_reconstruction_algorithm = the_global_data_dispatcher->get_reconstruction_algorithm();
}

bool track_finding_manager::find_tracks(vector<line_3d> cluster_lines[], vector<track>& the_tracks)
{
  int event_class;
  
  bool operation_success;
  
  double fit_goodness[max_number_of_tracks]; // this is a bad size of the array
  
  track upstream_track, downstream_track, total_track;
  
  double dummy_double;
  
  vector<pattern> pattern_options; // allow the possibility here, but assign_clusters_midstream should select a single one
  
  operation_success = classify_event(cluster_lines);

  if (!operation_success) return operation_success;
  
  event_class = the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity;
  
  if (!(event_class < 9000))
  {
    return false;
  }
  
  if (the_global_data_dispatcher->get_reconstruction_algorithm()->use_magnetic_field == false)
  {
    if (event_class == 1111 || event_class == 1112)
    {
      // event_class is 1111
      //
      // fit a straight line through all cluster lines; also fit two straight lines upstream and midstream plus downstream
      // compare which one is better; if the first one, then nothing happened; if the second one find the vertex;
      // if both are bad (by the chi square of the fit and the applied cut) either ignore or look for vertex on a plate (determined by option)
      // by fitting separately to and from that plate; the possibility exists that there can be a double vertex, which can be looked for too
      // (it can be target/plate or plate/plate case, the latter pretty rare), but perhaps ignore for now; this does not depend on alignment,
      // but alignment may affect errors and correlations
      //
      //
      // event_class is 1112
      //
      // fit a straight line through all the one-cluster plates, on the multiple cluster-planes add the one that is the closest (if both are really close use average?)
      // and fit again; do the same for midstream plus downstream; for upstream it can have one cluster on all (so it's done then), but it can have multiple clusters
      // on say one plate or zero clusters on say one plate; if multiple choose the one that produces a better vertex (smaller distance between the two fitted skew lines),
      // if zero determine the vertex by intersection of the plane formed by the other two lines and the midstream/downstream fit; this depends on the alignment of these
      // two strip lines, so will need to correct using the position on the third plate (such a procedure will be necessary in general); with multiple clusters upstream can
      // also use the one that will make the direction the closest to the beam direction
    
      // !!!! this needs to be reworked altogether !!!!
      
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
    
      fit_goodness[0] = upstream_track.fit_goodness;
    
      operation_success = find_midstream_plus_downstream_track_no_field(cluster_lines, downstream_track, event_class);
    
      fit_goodness[1] = downstream_track.fit_goodness;
    
      operation_success = find_total_track_no_field(cluster_lines, total_track, event_class);
    
      fit_goodness[2] = total_track.fit_goodness;
    
      XYZVector temp_vertex;
      find_closest_point(upstream_track.track_line, downstream_track.track_line, temp_vertex);
    
      //cout << "Fit goodness: " << fit_goodness[0] << " " << fit_goodness[1] << " " << fit_goodness[2] << endl;
      //cout << "Temp vertex: " << temp_vertex.X() << " " << temp_vertex.Y() << " " << temp_vertex.Z() << endl;
      
      /*
      if (event_class == 1111 && the_global_data_dispatcher->get_input_root_data_structure()->track_id->size() > 1)
      {
        cout << "from fit: " << upstream_track.track_line.direction << " " << downstream_track.track_line.direction << endl;
        cout << "-----------" << endl;
      }
      */
      
      //if (event_class == 1111 && the_global_data_dispatcher->get_input_root_data_structure()->track_id->size() > 1)
      /*
      if (event_class == 1111 && the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() > 1)
      {
        upstream_track.track_line.direction = upstream_track.track_line.direction.Unit();
        downstream_track.track_line.direction = downstream_track.track_line.direction.Unit();
      
        cout << "from fit: " << upstream_track.track_line.direction.Dot(downstream_track.track_line.direction) << endl;
      }
      */
      
      /*
      if (event_class == 1111)
      {
        cout << "--------------------" << endl;
        
        for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i++)
        {
          dummy_double = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->position.Z();
        
          cout << upstream_track.track_line.point.X() + (dummy_double - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.X() / upstream_track.track_line.direction.Z() << " "
               << upstream_track.track_line.point.Y() + (dummy_double - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.Y() / upstream_track.track_line.direction.Z() << " "
               << dummy_double << endl;
        }
        for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
        {
          dummy_double = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->position.Z();
          
          cout << downstream_track.track_line.point.X() + (dummy_double - downstream_track.track_line.point.Z()) * downstream_track.track_line.direction.X() / downstream_track.track_line.direction.Z() << " "
               << downstream_track.track_line.point.Y() + (dummy_double - downstream_track.track_line.point.Z()) * downstream_track.track_line.direction.Y() / downstream_track.track_line.direction.Z() << " "
               << dummy_double << endl;
        }
      }
      */
      
    }
    else if (event_class == 1117 || event_class == 2000)
    {
      // event_class is 1117
      //
      // look for one track midstream, but there are more than the allowed plates with multiple clusters (which
      // would qualify it for class 1112)
      //
      // event_class is 2000
      //
      // look for at least two tracks midstream
      
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
      
      if (operation_success)
      {
        operation_success = assign_clusters_midstream(cluster_lines, pattern_options, event_class, upstream_track);
      }
      
      if (operation_success)
      {
        operation_success = find_midstream_tracks_no_field(pattern_options, event_class, upstream_track);
      }
      
      /*
      // this is not needed in the no-field case, because there is no magnet
      if (operation_success)
      {
        operation_success = propagate_tracks_downstream(pattern_options);
      }
      */
      
      if (operation_success)
      {
        operation_success = assign_clusters_downstream_no_field(pattern_options, cluster_lines);
      }
cout << "aaaa" << endl;
      if (operation_success)
      {
        operation_success = do_global_fit_no_field(pattern_options, upstream_track);
      }
      
cout << "bbbb" << endl;

    }
  } // end if on magnetic field
  else
  {
    if (event_class == 1111 || event_class == 1112)
    {
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
      
      // assign clusters for 1112
      
      // propagate and fit or use Matej's code
    }
    else if (event_class == 1117 || event_class == 2000)
    {
      // find upstream track
      
      // assign clusters to tracks midstream
      
      // fit or calculate tracks midstream, assuming magnetic field has small effect
      
      // propagate through magnet and downstream
      
      // assign downstream clusters to tracks
      
      // propagate and fit or use Matej's code
      
      // find vertex
      
    }
  }
    
  return operation_success;
} 

bool track_finding_manager::classify_event(vector<line_3d> cluster_lines[])
{
  int temp_event_class = 0; // all event classes are positive integers
  
  bool is_9001, is_9002, is_9003, is_9004, is_9005;
  // not non-increasing (by new criteria), upstream too few plates hit, upstream too many plates with > 1 hit,
  // upstream too many hits on last plate (implying an early vertex)
  
  int number_of_plates, number_of_upstream_plates, number_of_midstream_plates, number_of_downstream_plates;
  int number_of_groups, number_of_upstream_groups, number_of_midstream_groups, number_of_downstream_groups;
  int number_of_upstream_XY_groups, number_of_midstream_XY_groups, number_of_downstream_XY_groups;
  int midstream_XY_groups[max_number_of_groups_in_tracking_region];
  int dummy_array[max_number_of_groups_in_tracking_region];
  int my_counter;
  
  int zero_cluster_counter_all_regions = 0;
  int zero_cluster_counter_for_region[max_number_of_tracking_regions];
  int single_cluster_counter_all_regions = 0;
  int single_cluster_counter_for_region[max_number_of_tracking_regions];
  int double_cluster_counter_all_regions = 0;
  int double_cluster_counter_for_region[max_number_of_tracking_regions];
  int multiple_cluster_counter_all_regions = 0;
  int multiple_cluster_counter_for_region[max_number_of_tracking_regions];
  
  int min_clusters;
  int dummy_int_1, dummy_int_2;
  
  int group_cluster_multiplicity[max_number_of_groups];
  
  // this may be a better indicator of where the vertex is
  int total_upstream_cluster_lines = 0;
  int total_midstream_cluster_lines = 0;
  int total_downstream_cluster_lines = 0;

  /*
  int number_of_xy_groups_midstream; // must be 2 or 3
  int xy_groups_midstream[max_number_of_groups_in_tracking_region];
  
  int number_of_dd_groups_midstream; // must be 0 or 1
  int dd_groups_midstream[max_number_of_groups_in_tracking_region];

  int number_of_diagonal_plates_midstream; // must be 1 or 2
  int diagonal_plates_midstream[2];
  */
  
  number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  number_of_midstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates();
  number_of_downstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_plates();
  
  number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  number_of_midstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups();
  number_of_downstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_groups();
  
  the_global_data_dispatcher->get_partial_geometry()->get_upstream_XY_groups(number_of_upstream_XY_groups, dummy_array);
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_XY_groups(number_of_midstream_XY_groups, midstream_XY_groups);
  the_global_data_dispatcher->get_partial_geometry()->get_downstream_XY_groups(number_of_downstream_XY_groups, dummy_array);
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    zero_cluster_counter_for_region[i] = 0;
    single_cluster_counter_for_region[i] = 0;
    double_cluster_counter_for_region[i] = 0;
    multiple_cluster_counter_for_region[i] = 0;
  }
  
  for (int i = 0; i < max_number_of_groups; i++)
  {
    group_cluster_multiplicity[i] = 0;
  }
  
  for (int i = 0; i < number_of_plates; i++)
  {
    if (i < number_of_upstream_plates)
    {
      total_upstream_cluster_lines += cluster_lines[i].size();
    }
    else if (i < number_of_upstream_plates + number_of_midstream_plates)
    {
      total_midstream_cluster_lines += cluster_lines[i].size();
    }
    else
    {
      total_downstream_cluster_lines += cluster_lines[i].size();
    }
    
    if (cluster_lines[i].size() == 0)
    {
      zero_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        zero_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        zero_cluster_counter_for_region[1]++;
      }
      else
      {
        zero_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() == 1)
    {
      single_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        single_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        single_cluster_counter_for_region[1]++;
      }
      else
      {
        single_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() == 2)
    {
      double_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        double_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        double_cluster_counter_for_region[1]++;
      }
      else
      {
        double_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() > 2)
    {
      multiple_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        multiple_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        multiple_cluster_counter_for_region[1]++;
      }
      else
      {
        multiple_cluster_counter_for_region[2]++;
      }
    }
  }
  
  // ?? this may need to change
  for (int i = 0; i < number_of_groups; i++)
  {
    min_clusters = 1000;
    
    dummy_int_1 = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->number_of_plates_in_group;
    
    for (int j = 0; j < dummy_int_1; j++)
    {
      dummy_int_2 = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->plate_index_list[j];
      
      if (the_global_data_dispatcher->get_partial_geometry()->get_plate(dummy_int_2)->intended_to_measure == "X" ||
          the_global_data_dispatcher->get_partial_geometry()->get_plate(dummy_int_2)->intended_to_measure == "Y")
      {
        if (cluster_lines[dummy_int_2].size() < min_clusters)
        {
          min_clusters = cluster_lines[dummy_int_2].size();
        }
      }
    }
    
    group_cluster_multiplicity[i] = min_clusters;
  }
  
  if (single_cluster_counter_all_regions == number_of_plates)
  {
    temp_event_class = 1111;
  }
  else if ((single_cluster_counter_for_region[0] + single_cluster_counter_for_region[1] + single_cluster_counter_for_region[2]) >=
           number_of_plates - the_global_data_dispatcher->get_reconstruction_options()->max_plates_with_multiple_clusters_for_single_track &&
           single_cluster_counter_for_region[1] >= number_of_midstream_plates - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plates_with_multiple_clusters_for_single_track &&
           single_cluster_counter_for_region[2] >= number_of_downstream_plates - the_global_data_dispatcher->get_reconstruction_options()->max_downstream_plates_with_multiple_clusters_for_single_track)
  {
    temp_event_class = 1112;
  }
  
  if (zero_cluster_counter_for_region[0] > the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_zero_clusters)
  {
    is_9002 = true;
  }
  else
  {
    is_9002 = false;
  }
  
  if (double_cluster_counter_for_region[0] + multiple_cluster_counter_for_region[0] > the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters)
  {
    is_9003 = true;
  }
  else
  {
    is_9003 = false;
  }
  
  if (cluster_lines[number_of_upstream_plates - 1].size() > 2) // make it a cut; maybe 3 is better
  {
    is_9004 = true;
  }
  else
  {
    is_9004 = false;
  }
  
  if (total_upstream_cluster_lines > 6) // make it a cut
  {
    is_9005 = true;
  }
  else
  {
    is_9005 = false;
  }
  
  if ( (((double)total_midstream_cluster_lines) / ((double)(number_of_midstream_plates))) /
       (((double)total_upstream_cluster_lines) / ((double)(number_of_upstream_plates))) >= 1.6 &&
         cluster_lines[number_of_upstream_plates].size() > 1)
    // the number 1.6 needs to be a cut
    // 1.6 includes 1 1 1 2  2 2 2 2 2
  {
    is_9001 = false;
  }
  else
  {
    if (temp_event_class < 1000)
    {
      is_9001 = true;
    }
    else
    {
      is_9001 = false;
    }
  }
  
  my_counter = 0;
  
  for (int i = 0; i < number_of_midstream_XY_groups; i++)
  {
    if (group_cluster_multiplicity[midstream_XY_groups[i]] >= 2)
    {
      my_counter++;
    }
  }
  
  if (my_counter < 2 && temp_event_class < 1000) // make 2 a cut
  {
    temp_event_class = 1117;
  }
  
  /*
  for (int i = 0; i < number_of_plates; i++)
  {
    cout << cluster_lines[i].size() << " ";
  }
  cout << endl;
  */
   
  /*
  if (is_9001) cout << 9001 << " ";
  if (is_9002) cout << 9002 << " ";
  if (is_9003) cout << 9003 << " ";
  if (is_9004) cout << 9004 << " ";
  if (is_9005) cout << 9005 << " ";
  if (temp_event_class == 1111) cout << 1111 << " " ;
  if (temp_event_class == 1112) cout << 1112 << " " ;
  if (temp_event_class == 1117) cout << 1117 << " " ;
  cout << endl;
  */
  
  if (temp_event_class != 1111 && temp_event_class != 1112 && temp_event_class != 1117) temp_event_class = 2000;
    
  if (is_9001 || is_9002 || is_9003 || is_9004 || is_9005) temp_event_class = 9000;
  
  // cout << "Event class: " << temp_event_class << endl;
  
  
  
  
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity = temp_event_class;
  
  
  
  return true;
}

bool track_finding_manager::find_upstream_track(vector<line_3d> cluster_lines[], track& upstream_track, int event_class)
{
  bool operation_success;
  
  int plate_index_multiple_clusters;
  int plate_index_zero_clusters;
  
  track temp_track, temp_best_track;
  double best_direction, temp_direction;
  XYZVector z_direction;
  
  string zero_cluster_plate_type;
  line_3d temp_cluster_line;
  
  // general enough ??
  z_direction.SetCoordinates(0., 0., 1.);
  
  int count_zero_cluster_plates = 0;
  int count_multiple_cluster_plates = 0;
  
  int number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  
  upstream_track.cluster_lines.clear();
  upstream_track.plate_index_of_cluster_line.clear();
  
  for (int i = 0; i < number_of_upstream_plates; i++)
  {
    if (cluster_lines[i].size() == 1)
    {
      upstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
      upstream_track.plate_index_of_cluster_line.push_back(i);
    }
    else if (cluster_lines[i].size() == 0)
    {
      count_zero_cluster_plates++;
      plate_index_zero_clusters = i;
    }
    else
    {
      count_multiple_cluster_plates++;
      plate_index_multiple_clusters = i;
    }
  }
  
  // check cuts; !! do properly later
  if (count_zero_cluster_plates + count_multiple_cluster_plates > 1)
  {
    cout << "Problem: too many upstream plates with cluster number different than 1." << endl;
    return false;
  }
  else if (count_zero_cluster_plates + count_multiple_cluster_plates == 1 && event_class == 1111)
  {
    cout << "Problem: event class 1111 with clusters != 1 on an upstream plate." << endl;
    return false;
  }
  
  upstream_track.assumed_start_group = 0;
  upstream_track.assumed_end_group = number_of_upstream_groups - 1;
  
  if (count_zero_cluster_plates + count_multiple_cluster_plates == 0)
  {
    operation_success = fit_track_no_field(upstream_track);
    // calculate the position on target and error (both at z = 0)
    return operation_success;
  }
  else if (count_multiple_cluster_plates == 1)
  {
    best_direction = 10.;
    
    for (int i = 0; i < cluster_lines[plate_index_multiple_clusters].size(); i++)
    {
      temp_track.cluster_lines.clear();
      temp_track.plate_index_of_cluster_line.clear();
      
      temp_track = upstream_track;
      temp_track.cluster_lines.push_back(cluster_lines[plate_index_multiple_clusters].at(i));
      temp_track.plate_index_of_cluster_line.push_back(plate_index_multiple_clusters);
      
      temp_track.sort_cluster_lines_by_plate_index();
      operation_success = fit_track_no_field(temp_track);
      // calculate the position on target and error (both at z = 0)
      temp_direction = fabs(sqrt(temp_track.track_line.direction.Unit().Cross(z_direction).Mag2()));
      
      if (temp_direction < best_direction)
      {
        best_direction = temp_direction;
        temp_best_track = temp_track;
      }
    }
    
    // do with cut later !!
    if (best_direction < 0.01)
    {
      upstream_track = temp_best_track;
      return true;
    }
    else
    {
      return false;
    }
  }
  else if (count_zero_cluster_plates == 1)
  {
    zero_cluster_plate_type = the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index_zero_clusters)->intended_to_measure;
    
    for (int i = 0; i < number_of_upstream_plates; i++)
    {
      if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == zero_cluster_plate_type)
      {
        for (int j = 0; j < upstream_track.cluster_lines.size(); j++)
        {
          if (upstream_track.plate_index_of_cluster_line.at(j) == i)
          {
            temp_cluster_line.point = upstream_track.cluster_lines.at(j).point;
            temp_cluster_line.point.SetZ(the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index_zero_clusters)->position.Z());
            temp_cluster_line.direction = upstream_track.cluster_lines.at(j).direction;
          }
        }
      }
    }
    
    upstream_track.cluster_lines.push_back(temp_cluster_line);
    upstream_track.plate_index_of_cluster_line.push_back(plate_index_zero_clusters);
    upstream_track.sort_cluster_lines_by_plate_index();
    
    operation_success = fit_track_no_field(upstream_track);
    // calculate the position on target and error (both at z = 0)
    
    return operation_success;
  }
  else
  {
    cout << "A track_finding_manager::find_upstream_track exception ??. This shouldn't happen...";
    return false;
  }
}

bool track_finding_manager::find_midstream_plus_downstream_track_no_field(std::vector<line_3d> cluster_lines[], track& middownstream_track, int event_class)
{
  bool operation_success;
  
  int number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  int number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  int number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  vector<int> unused_plates_list;
  int closest_cluster;
  double closest_cluster_distance, dummy_distance;
  
  middownstream_track.cluster_lines.clear();
  middownstream_track.plate_index_of_cluster_line.clear();
  
  middownstream_track.assumed_start_group = number_of_upstream_groups;
  middownstream_track.assumed_end_group = number_of_groups - 1;
  
  if (event_class == 1111)
  {
    for (int i = number_of_upstream_plates; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
        middownstream_track.plate_index_of_cluster_line.push_back(i);
      }
      else
      {
        cout << "Problem: event class 1111 with clusters != 1 on a midstream/downstream plate." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    operation_success = fit_track_no_field(middownstream_track);
    
    return operation_success;
  }
  else if (event_class == 1112)
  {
    middownstream_track.cluster_lines.clear();
    middownstream_track.plate_index_of_cluster_line.clear();
    unused_plates_list.clear();
    
    int primary_direction_1_plate_counter = 0;
    int primary_direction_2_plate_counter = 0;
    
    for (int i = number_of_upstream_plates; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
        middownstream_track.plate_index_of_cluster_line.push_back(i);
        
        if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name())
        {
          primary_direction_1_plate_counter++;
        }
        else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name())
        {
          primary_direction_2_plate_counter++;
        }
      }
      else if (cluster_lines[i].size() > 1)
      {
        unused_plates_list.push_back(i);
      }
    }
    
    if (primary_direction_1_plate_counter < 2 || primary_direction_2_plate_counter < 2)
    {
      cout << "Track Finding Manager Error ???? See code. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      operation_success = fit_track_no_field(middownstream_track);
    }
    
    for (int i = 0; i < unused_plates_list.size(); i++)
    {
      closest_cluster = -1;
      closest_cluster_distance = 1000.;
      
      if (cluster_lines[unused_plates_list.at(i)].size() > 1)
      {
        for (int j = 0; j < cluster_lines[unused_plates_list.at(i)].size(); j++)
        {
          if (find_line_to_line_distance(cluster_lines[unused_plates_list.at(i)].at(j), middownstream_track.track_line, dummy_distance))
          {
            if (dummy_distance < closest_cluster_distance)
            {
              closest_cluster = j;
              closest_cluster_distance = dummy_distance;
            }
          }
        }
      }
      
      if (closest_cluster_distance < 0.01) // !! use cut later
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[unused_plates_list.at(i)].at(closest_cluster));
        middownstream_track.plate_index_of_cluster_line.push_back(unused_plates_list.at(i));
      }
    }
    
    middownstream_track.sort_cluster_lines_by_plate_index();
    
    operation_success = fit_track_no_field(middownstream_track);
    
    return operation_success;
  }
}

bool track_finding_manager::find_total_track_no_field(std::vector<line_3d> cluster_lines[], track& total_track, int event_class)
{
  bool operation_success;
  
  vector<int> unused_plates_list;
  int closest_cluster;
  double closest_cluster_distance, dummy_distance;
  
  int number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  int number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  total_track.cluster_lines.clear();
  total_track.plate_index_of_cluster_line.clear();
  
  total_track.assumed_start_group = 0;
  total_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups() - 1;
  
  if (event_class == 1111)
  {
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        total_track.cluster_lines.push_back(cluster_lines[i].at(0));
        total_track.plate_index_of_cluster_line.push_back(i);
      }
      else
      {
        cout << "Problem: event class 1111 with clusters != 1 on an upstream/midstream/downstream plate." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    operation_success = fit_track_no_field(total_track);
    return operation_success;
  }
  else if (event_class == 1112)
  {
    total_track.cluster_lines.clear();
    total_track.plate_index_of_cluster_line.clear();
    unused_plates_list.clear();
    
    // used to make sure there are at least 2 plates in each primary direction for the initial fit
    int primary_direction_1_plate_counter = 0;
    int primary_direction_2_plate_counter = 0;
    
    for (int i = 0; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        total_track.cluster_lines.push_back(cluster_lines[i].at(0));
        total_track.plate_index_of_cluster_line.push_back(i);
        
        if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name())
        {
          primary_direction_1_plate_counter++;
        }
        else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name())
        {
          primary_direction_2_plate_counter++;
        }
      }
      else if (cluster_lines[i].size() > 1)
      {
        unused_plates_list.push_back(i);
      }
    }
    
    if (primary_direction_1_plate_counter < 2 || primary_direction_2_plate_counter < 2)
    {
      cout << "Track Finding Manager Error ???? See code. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      operation_success = fit_track_no_field(total_track);
    }
    
    for (int i = 0; i < unused_plates_list.size(); i++)
    {
      closest_cluster = -1;
      closest_cluster_distance = 1000.;
      
      if (cluster_lines[unused_plates_list.at(i)].size() > 1)
      {
        for (int j = 0; j < cluster_lines[unused_plates_list.at(i)].size(); j++)
        {
          if (find_line_to_line_distance(cluster_lines[unused_plates_list.at(i)].at(j), total_track.track_line, dummy_distance))
          {
            if (dummy_distance < closest_cluster_distance)
            {
              closest_cluster = j;
              closest_cluster_distance = dummy_distance;
            }
          }
        }
      }
      
      if (closest_cluster_distance < 0.01) // !! use cut later
      {
        total_track.cluster_lines.push_back(cluster_lines[unused_plates_list.at(i)].at(closest_cluster));
        total_track.plate_index_of_cluster_line.push_back(unused_plates_list.at(i));
      }
    }
    
    total_track.sort_cluster_lines_by_plate_index();
    operation_success = fit_track_no_field(total_track);
    
    return operation_success;
  }
}

bool track_finding_manager::fit_track_no_field(track& a_track)
{
  bool fitting_success = true;
  
  double average_x = 0., average_y = 0., average_z = 0.;
  double dx_dz, dy_dz;
  
  string primary_direciton_names[2];
  int primary_direction_plates[2][max_number_of_plates];
  int number_of_primary_direction_plates[2] = {0, 0};
  
  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < max_number_of_plates; j++)
    {
      primary_direction_plates[i][j] = -1;
    }
  }
  
  primary_direciton_names[0] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name();
  primary_direciton_names[1] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name();
  
  for (int i = 0; i < a_track.cluster_lines.size(); i++)
  {
    if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->intended_to_measure == primary_direciton_names[0])
    {
      primary_direction_plates[0][number_of_primary_direction_plates[0]] = i;
      number_of_primary_direction_plates[0]++;
    }
    else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->intended_to_measure == primary_direciton_names[1])
    {
      primary_direction_plates[1][number_of_primary_direction_plates[1]] = i;
      number_of_primary_direction_plates[1]++;
    }
  }
  
  if (number_of_primary_direction_plates[0] < 2 || number_of_primary_direction_plates[1] < 2)
  {
    cout << "Too few plates in a given direction to fit a track..." << endl;
    return false;
  }
  
  // also only implemented for X and Y as primary
  for (int i = 0; i < number_of_primary_direction_plates[0]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.z();
    
    if (primary_direciton_names[0] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.x();
    }
    else if (primary_direciton_names[0] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[0] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[0];
  }
  else if (primary_direciton_names[0] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[0];
  }
  
  for (int i = 0; i < number_of_primary_direction_plates[1]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.z();
    
    if (primary_direciton_names[1] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.x();
    }
    else if (primary_direciton_names[1] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[1] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[1];
  }
  else if (primary_direciton_names[1] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[1];
  }
  
  average_z /= (double)(number_of_primary_direction_plates[0] + number_of_primary_direction_plates[1]);
  
  // this is a particular case; work out the general case later
  //
  if (primary_direciton_names[0] == "X" && primary_direciton_names[1] == "Y")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
  }
  else if (primary_direciton_names[0] == "Y" && primary_direciton_names[1] == "X")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
  }
  else
  {
    cout << "This configuration of primary directions has not been implemented yet..." << endl;
    return false;
  }
  
  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();
  
  initial_values_of_parameters.push_back(average_x);
  initial_values_of_parameters.push_back(average_y);
  initial_values_of_parameters.push_back(dx_dz);
  initial_values_of_parameters.push_back(dy_dz);
  
  // need better estimates for these
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.01);
  errors_in_initial_values_of_parameters.push_back(0.01);
  
  my_straight_line_fcn the_fcn(a_track.cluster_lines, average_z);
  
  VariableMetricMinimizer the_minimizer;
  
  FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
  a_track.track_line.error_point.SetCoordinates(min.UserState().Error(0), min.UserState().Error(1), 0.);
  // dangerous, may lead to division by 0 !!
  a_track.track_line.error_direction.SetCoordinates(min.UserState().Error(2) / min.UserState().Value(2), min.UserState().Error(3) / min.UserState().Value(3), 0.);
  a_track.track_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), average_z);
  a_track.track_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
  a_track.track_line.direction = a_track.track_line.direction.Unit();
  
  a_track.fit_goodness = min.Fval();
  
  return fitting_success;
}

bool track_finding_manager::assign_clusters_midstream(vector<line_3d> cluster_lines[], vector<pattern>& pattern_options, int event_class, track upstream_track)
{
  bool operation_success = true;
  
  pattern temp_pattern;
  vector<pattern> x_pattern_options, y_pattern_options;
  track temp_track;
  
  // these are intended to measure
  int number_of_x_plates, number_of_y_plates, number_of_d1_plates, number_of_d2_plates, number_of_d_plates;
  ssd_plate x_plates[max_number_of_plates_in_tracking_region], y_plates[max_number_of_plates_in_tracking_region], d1_plates[max_number_of_plates_in_tracking_region], d2_plates[max_number_of_plates_in_tracking_region], d_plates[max_number_of_plates_in_tracking_region];
  int x_plate_indices[max_number_of_plates_in_tracking_region],  y_plate_indices[max_number_of_plates_in_tracking_region], d1_plate_indices[max_number_of_plates_in_tracking_region], d2_plate_indices[max_number_of_plates_in_tracking_region], d_plate_indices[max_number_of_plates_in_tracking_region];

  double midstream_X_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_X_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  
  int u_plate, d_plate; // up- and downstream of the two plates
  int n_u_clusters, n_d_clusters; // number of u- and d-plate clusters
  
  double scale_for_2;
  double u_scale, d_scale;
  double cut_min_tolerance = 0., cut_max_tolerance = 0.;

  double start_z, stop_z, step_z, current_z;
  int corresponding_u_points[200], best_corresponding_u_points[200]; // this should be max_clusters_per_plate
  double distances[200]; // ditto
  XYZVector target_point, best_target_point;
  double projection_coordinates[200]; // ditto
  double total_distance, best_total_distance;
  
  line_3d temp_line;
  
  double option_best_distance[2];
  int one_track_options[6];
  int two_track_options[2][6];
  double min_distance_on_diagonal_plates[2];
  double best_dummy_distance;
  double dummy_distance;
  double multi_track_best_distance[200][200]; // use cut later
  int multi_track_closest_d_strip[200][200][2]; // ditto
  int temp_closest_d_strip[2];
  double total_diagonal_distance;
  vector<multi_track_selections> the_multi_track_selections;

  line_3d one_track_line;
  line_3d two_track_line[2];
  line_3d multi_track_line[200][200]; // use cut later
  
  
  
  
  
  the_global_data_dispatcher->get_partial_geometry()->get_types_of_midstream_plates(number_of_x_plates, x_plates, x_plate_indices, number_of_y_plates, y_plates, y_plate_indices, number_of_d1_plates, d1_plates, d1_plate_indices, number_of_d2_plates, d2_plates, d2_plate_indices, number_of_d_plates, d_plates, d_plate_indices);
  
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_X_scales(midstream_X_scales_max, midstream_X_scales_min);
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_Y_scales(midstream_Y_scales_max, midstream_Y_scales_min);
  
  // do for both directions; Y goes first because it is not affected by the magnetic field
  // and is therefore more important
  if (number_of_y_plates == 2) // this part of the algorithm will not work if there are more than 2 X or 2 Y plates
  {
    u_plate = y_plate_indices[0];
    d_plate = y_plate_indices[1];
    
    n_u_clusters = cluster_lines[y_plate_indices[0]].size();
    n_d_clusters = cluster_lines[y_plate_indices[1]].size();
cout << "Y n_u, n_d " << n_u_clusters << " " << n_d_clusters << endl;
    if (n_u_clusters == 1 && n_d_clusters == 1)
    {
      y_pattern_options.clear();
      temp_pattern.reset();

      temp_track.cluster_lines.clear();
      temp_track.plate_index_of_cluster_line.clear();
      temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
      temp_track.plate_index_of_cluster_line.push_back(u_plate);
      temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
      temp_track.plate_index_of_cluster_line.push_back(d_plate);
      
      temp_pattern.y_tracks.push_back(temp_track);
      temp_pattern.number_of_y_lines = 1;
      y_pattern_options.push_back(temp_pattern);
    }
    else if (n_u_clusters > 1 && n_d_clusters == 1)
    {
      y_pattern_options.clear();
      
      for (int i = 0; i < n_u_clusters; i++)
      {
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.y_tracks.push_back(temp_track);
        temp_pattern.number_of_y_lines = 1;
        y_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters == 1 && n_d_clusters > 1)
    {
      y_pattern_options.clear();
      
      for (int i = 0; i < n_d_clusters; i++)
      {
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.y_tracks.push_back(temp_track);
        temp_pattern.number_of_y_lines = 1;
        y_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters == 2 && n_d_clusters == 2)
    {
      scale_for_2 = fabs((cluster_lines[d_plate].at(1).point.Y() - cluster_lines[d_plate].at(0).point.Y()) / (cluster_lines[u_plate].at(1).point.Y() - cluster_lines[u_plate].at(0).point.Y()));
      
      // if more than 2 Y plates the scales are more than 1 number, so need to modify code
      if (scale_for_2 > midstream_Y_scales_min[0] - cut_min_tolerance && scale_for_2 < midstream_Y_scales_max[0] + cut_max_tolerance)
      {
        y_pattern_options.clear();
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.y_tracks.push_back(temp_track);

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.y_tracks.push_back(temp_track);
        temp_pattern.number_of_y_lines = 2;
        y_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters > 2 && n_d_clusters == 2)
    {
      d_scale = fabs(cluster_lines[d_plate].at(1).point.Y() - cluster_lines[d_plate].at(0).point.Y());
      
      y_pattern_options.clear();

      for (int i = 0; i < n_u_clusters - 1; i++)
      {
        for (int j = n_u_clusters - 1; j > i; j--)
        {
          u_scale = fabs(cluster_lines[u_plate].at(j).point.Y() - cluster_lines[u_plate].at(i).point.Y());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_Y_scales_min[0] && d_scale / u_scale < midstream_Y_scales_max[0])
          {
            temp_pattern.reset();

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.y_tracks.push_back(temp_track);

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(j));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.y_tracks.push_back(temp_track);
            temp_pattern.number_of_y_lines = 2;
            y_pattern_options.push_back(temp_pattern);
          }
        }
      }
      /*
      if (y_tracks.size() == 0)
      {
        success_y_scale = false;
      }
      else
      {
        success_y_scale = true;
      }
      */
    }
    else if (n_u_clusters == 2 && n_d_clusters > 2)
    {
      u_scale = fabs(cluster_lines[u_plate].at(1).point.Y() - cluster_lines[u_plate].at(0).point.Y());
      
      y_pattern_options.clear();

      for (int i = 0; i < n_d_clusters - 1; i++)
      {
        for (int j = n_d_clusters - 1; j > i; j--)
        {
          d_scale = fabs(cluster_lines[d_plate].at(j).point.Y() - cluster_lines[d_plate].at(i).point.Y());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_Y_scales_min[0] && d_scale / u_scale < midstream_Y_scales_max[0])
          {
            temp_pattern.reset();

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.y_tracks.push_back(temp_track);

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(j));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.y_tracks.push_back(temp_track);
            temp_pattern.number_of_y_lines = 2;
            y_pattern_options.push_back(temp_pattern);
          }
        }
      }
      /*
      if (y_tracks.size() == 0)
      {
        success_y_scale = false;
      }
      else
      {
        success_y_scale = true;
      }
      */
    }
    else if (n_u_clusters > 2 && n_d_clusters > 2)
    {
      const int number_of_divisions = 100;
      
      start_z = the_global_data_dispatcher->get_partial_geometry()->get_target_min_z();
      stop_z = the_global_data_dispatcher->get_partial_geometry()->get_target_max_z();
      step_z = (stop_z - start_z) / ((double)(number_of_divisions - 1));
      
      best_total_distance = 1000.;
      
      y_pattern_options.clear();

      //cout << "------------------" << endl;
      for (int i = 0; i < number_of_divisions; i++)
      {
        current_z = start_z + step_z * ((double) i);
        
        target_point = find_point_on_line_by_z(upstream_track.track_line, current_z);
        //cout << target_point << endl;
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          corresponding_u_points[j] = -1;
          distances[j] = 100.;
          
          projection_coordinates[j] = target_point.Y() + (cluster_lines[d_plate].at(j).point.Y() - target_point.Y()) * (y_plates[0].position.Z() - target_point.Z()) / (y_plates[1].position.Z() - target_point.Z());
        }
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          for (int k = 0; k < n_u_clusters; k++)
          {
            //cout << "-->" << fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.Y()) << endl;
            if (fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.Y()) < distances[j])
            {
              distances[j] = fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.Y());
              corresponding_u_points[j] = k;
            }
          }
          //cout << "->" << distances[j] << endl;
        }
        
        total_distance = 0.;
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          total_distance += distances[j];
        }
        //cout << ">" << total_distance << endl;
        
        if (total_distance < best_total_distance)
        {
          best_total_distance = total_distance;
          best_target_point = target_point;
          
          for (int j = 0; j < n_d_clusters; j++)
          {
            best_corresponding_u_points[j] = corresponding_u_points[j];
          }
        }
      }
      
      the_global_data_dispatcher->get_event_reco_output_data_structure()->guessed_vertex = best_target_point;
      
       //cout << best_target_point << " " << best_total_distance << " ";
      for (int j = 0; j < n_d_clusters; j++)
      {
        // cout << best_corresponding_u_points[j] << " ";
      }
      // cout << endl;
      
      temp_pattern.reset();

      for (int i = 0; i < n_d_clusters; i++)
      {
        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(corresponding_u_points[i]));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.y_tracks.push_back(temp_track);
      }
      
      y_pattern_options.push_back(temp_pattern);
    }
    else
    {
      cout << "Not enough clusters on a plate with presumed at least 1 cluster. " << endl;
      // success_y_scale = false;
    }
  }
  else // more than 2 x or y planes midstream can be done later; not used in 2021 run
  {
    cout << "Midstream Y plates different than 2. Can't do cluster assignment. Quitting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  // now do the X plates
  if (number_of_x_plates == 2)
  {
    u_plate = x_plate_indices[0];
    d_plate = x_plate_indices[1];

    n_u_clusters = cluster_lines[x_plate_indices[0]].size();
    n_d_clusters = cluster_lines[x_plate_indices[1]].size();
cout << "X n_u, n_d " << n_u_clusters << " " << n_d_clusters << endl;
    if (n_u_clusters == 1 && n_d_clusters == 1)
    {
      x_pattern_options.clear();
      temp_pattern.reset();

      temp_track.cluster_lines.clear();
      temp_track.plate_index_of_cluster_line.clear();
      temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
      temp_track.plate_index_of_cluster_line.push_back(u_plate);
      temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
      temp_track.plate_index_of_cluster_line.push_back(d_plate);
      
      temp_pattern.x_tracks.push_back(temp_track);
      temp_pattern.number_of_x_lines = 1;
      x_pattern_options.push_back(temp_pattern);
    }
    else if (n_u_clusters > 1 && n_d_clusters == 1)
    {
      x_pattern_options.clear();
      
      for (int i = 0; i < n_u_clusters; i++)
      {
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.x_tracks.push_back(temp_track);
        temp_pattern.number_of_x_lines = 1;
        x_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters == 1 && n_d_clusters > 1)
    {
      x_pattern_options.clear();
      
      for (int i = 0; i < n_d_clusters; i++)
      {
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.x_tracks.push_back(temp_track);
        temp_pattern.number_of_x_lines = 1;
        x_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters == 2 && n_d_clusters == 2)
    {
      //cout << "Y22" << endl;
      
      scale_for_2 = fabs((cluster_lines[d_plate].at(1).point.X() - cluster_lines[d_plate].at(0).point.X()) / (cluster_lines[u_plate].at(1).point.X() - cluster_lines[u_plate].at(0).point.X()));
      
      //cout << scale_for_2 << endl;
      
      // if more than 2 Y plates the scales are more than 1 number, so need to modify code
      if (scale_for_2 > midstream_X_scales_min[0] - cut_min_tolerance && scale_for_2 < midstream_X_scales_max[0] + cut_max_tolerance)
      {
        x_pattern_options.clear();
        temp_pattern.reset();

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.x_tracks.push_back(temp_track);

        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        temp_pattern.x_tracks.push_back(temp_track);
        temp_pattern.number_of_x_lines = 2;
        x_pattern_options.push_back(temp_pattern);
      }
    }
    else if (n_u_clusters > 2 && n_d_clusters == 2)
    {
      d_scale = fabs(cluster_lines[d_plate].at(1).point.X() - cluster_lines[d_plate].at(0).point.X());
      
      x_pattern_options.clear();

      for (int i = 0; i < n_u_clusters - 1; i++)
      {
        for (int j = n_u_clusters - 1; j > i; j--)
        {
          u_scale = fabs(cluster_lines[u_plate].at(j).point.X() - cluster_lines[u_plate].at(i).point.X());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_X_scales_min[0] && d_scale / u_scale < midstream_X_scales_max[0])
          {
            temp_pattern.reset();

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.x_tracks.push_back(temp_track);

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(j));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);

            temp_pattern.x_tracks.push_back(temp_track);
            temp_pattern.number_of_x_lines = 2;
            x_pattern_options.push_back(temp_pattern);          }
        }
      }
      /*
      if (x_tracks.size() == 0)
      {
        success_x_scale = false;
      }
      else
      {
        success_x_scale = true;
      }
      */
    }
    else if (n_u_clusters == 2 && n_d_clusters > 2)
    {
      u_scale = fabs(cluster_lines[u_plate].at(1).point.X() - cluster_lines[u_plate].at(0).point.X());
      
      x_pattern_options.clear();

      for (int i = 0; i < n_d_clusters - 1; i++)
      {
        for (int j = n_d_clusters - 1; j > i; j--)
        {
          d_scale = fabs(cluster_lines[d_plate].at(j).point.X() - cluster_lines[d_plate].at(i).point.X());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_X_scales_min[0] && d_scale / u_scale < midstream_X_scales_max[0])
          {
            temp_pattern.reset();

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.x_tracks.push_back(temp_track);

            temp_track.cluster_lines.clear();
            temp_track.plate_index_of_cluster_line.clear();
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
            temp_track.plate_index_of_cluster_line.push_back(u_plate);
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(j));
            temp_track.plate_index_of_cluster_line.push_back(d_plate);
            
            temp_pattern.x_tracks.push_back(temp_track);
            temp_pattern.number_of_x_lines = 2;
            x_pattern_options.push_back(temp_pattern);
          }
        }
      }
      /*
      if (x_tracks.size() == 0)
      {
        success_x_scale = false;
      }
      else
      {
        success_x_scale = true;
      }
      */
    }
    else if (n_u_clusters > 2 && n_d_clusters > 2)
    {
      const int number_of_divisions = 100;
      
      start_z = the_global_data_dispatcher->get_partial_geometry()->get_target_min_z();
      stop_z = the_global_data_dispatcher->get_partial_geometry()->get_target_max_z();
      step_z = (stop_z - start_z) / ((double)(number_of_divisions - 1));
      
      best_total_distance = 1000.;
      
      x_pattern_options.clear();

      //cout << "------------------" << endl;
      for (int i = 0; i < number_of_divisions; i++)
      {
        current_z = start_z + step_z * ((double) i);
        
        target_point = find_point_on_line_by_z(upstream_track.track_line, current_z);
        //cout << target_point << endl;
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          corresponding_u_points[j] = -1;
          distances[j] = 100.;
          
          projection_coordinates[j] = target_point.X() + (cluster_lines[d_plate].at(j).point.X() - target_point.X()) * (x_plates[0].position.Z() - target_point.Z()) / (x_plates[1].position.Z() - target_point.Z());
        }
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          for (int k = 0; k < n_u_clusters; k++)
          {
            //cout << "-->" << fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.X()) << endl;
            if (fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.X()) < distances[j])
            {
              distances[j] = fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.X());
              corresponding_u_points[j] = k;
            }
          }
          //cout << "->" << distances[j] << endl;
        }
        
        total_distance = 0.;
        
        for (int j = 0; j < n_d_clusters; j++)
        {
          total_distance += distances[j];
        }
        //cout << ">" << total_distance << endl;
        
        if (total_distance < best_total_distance)
        {
          best_total_distance = total_distance;
          best_target_point = target_point;
          
          for (int j = 0; j < n_d_clusters; j++)
          {
            best_corresponding_u_points[j] = corresponding_u_points[j];
          }
        }
      }
      
      //cout << best_target_point << " " << best_total_distance << " ";
      for (int j = 0; j < n_d_clusters; j++)
      {
        //cout << best_corresponding_u_points[j] << " ";
      }
      //cout << endl;
      
      temp_pattern.reset();

      for (int i = 0; i < n_d_clusters; i++)
      {
        temp_track.cluster_lines.clear();
        temp_track.plate_index_of_cluster_line.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(corresponding_u_points[i]));
        temp_track.plate_index_of_cluster_line.push_back(u_plate);
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
        temp_track.plate_index_of_cluster_line.push_back(d_plate);
        
        if (n_u_clusters < n_d_clusters)
        {
          temp_pattern.number_of_x_lines = n_u_clusters;
        }
        else
        {
          temp_pattern.number_of_x_lines = n_d_clusters;
        }
          
        temp_pattern.x_tracks.push_back(temp_track);
      }
      
      x_pattern_options.push_back(temp_pattern);
    }
    else
    {
      cout << "Not enough clusters on a plate with presumed at least 1 cluster. " << endl;
    }
  }
  else // more than 2 x or y planes midstream can be done later; not used in 2021 run
  {
    cout << "Midstream X plates different than 2. Can't do cluster assignment. Quitting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  // now combine the X tracks and the Y tracks in all possible combinations
  //
  pattern_options.clear();

  option_best_distance[0] = 100.;
  option_best_distance[1] = 100.;
  
  one_track_options[0] = -1;
  one_track_options[1] = -1;
  one_track_options[2] = -1;
  one_track_options[3] = -1;
  one_track_options[4] = -1;
  one_track_options[5] = -1;
  
  two_track_options[0][0] = -1;
  two_track_options[0][1] = -1;
  two_track_options[0][2] = -1;
  two_track_options[0][3] = -1;
  two_track_options[0][4] = -1;
  two_track_options[0][5] = -1;
  two_track_options[1][0] = -1;
  two_track_options[1][1] = -1;
  two_track_options[1][2] = -1;
  two_track_options[1][3] = -1;
  two_track_options[1][4] = -1;
  two_track_options[1][5] = -1;
  
  for (int i = 0; i < 200; i++) // use cuts later
  {
    for (int j = 0; j < 200; j++)
    {
      multi_track_best_distance[i][j] = 100.;
      multi_track_closest_d_strip[i][j][0] = -1;
      multi_track_closest_d_strip[i][j][1] = -1;
    }
  }
  
  temp_closest_d_strip[0] = -1;
  temp_closest_d_strip[1] = -1;
  
  
cout << "----------------------" << endl;
cout << "x and y options: " << x_pattern_options.size() << " " << y_pattern_options.size() << endl;
  for (int i = 0; i < x_pattern_options.size(); i++)
  {
    for (int j = 0; j < y_pattern_options.size(); j++)
    {
cout << "x and y tracks: " << x_pattern_options.at(i).x_tracks.size() << " " << y_pattern_options.at(j).y_tracks.size() << endl;

      option_best_distance[0] = 100.;
      option_best_distance[1] = 100.;
      
      if (x_pattern_options.at(i).x_tracks.size() == 1 || y_pattern_options.at(j).y_tracks.size() == 1)
      {
        for (int m = 0; m < x_pattern_options.at(i).x_tracks.size(); m++)
        {
          for (int n = 0; n < y_pattern_options.at(j).y_tracks.size(); n++)
          {
            
            // !!!! this is assuming there are 2 X plates and 2 Y plates;
            // if there are more the function needs to be modified
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(m), y_pattern_options.at(j).y_tracks.at(n), temp_line);
            
            min_distance_on_diagonal_plates[0] = 0.;
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[k] = kk;
                }
              }
              
              min_distance_on_diagonal_plates[0] += best_dummy_distance;
            }
            
            if (min_distance_on_diagonal_plates[0] < option_best_distance[0])
            {
              option_best_distance[0] = min_distance_on_diagonal_plates[0];
              
              one_track_options[0] = i;
              one_track_options[1] = j;
              one_track_options[2] = m;
              one_track_options[3] = n;
              one_track_options[4] = temp_closest_d_strip[0];
              one_track_options[5] = temp_closest_d_strip[1];
              one_track_line = temp_line;
            }
          }
        }
        
        if (option_best_distance[0] < 0.2) // !! make cut later
        {
          temp_pattern.reset();
          
          temp_pattern.pattern_quality[0] = option_best_distance[0];
          
          temp_track.cluster_lines.clear();
          temp_track.plate_index_of_cluster_line.clear();
          
          for (int k = 0; k < x_pattern_options.at(one_track_options[0]).x_tracks.at(one_track_options[2]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(x_pattern_options.at(one_track_options[0]).x_tracks.at(one_track_options[2]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(x_pattern_options.at(one_track_options[0]).x_tracks.at(one_track_options[2]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[0]].at(one_track_options[4]));
          temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[0]);
          if (number_of_d_plates > 1)
          {
            temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[1]].at(one_track_options[5]));
            temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[1]);
          }
          
          for (int k = 0; k < y_pattern_options.at(one_track_options[1]).y_tracks.at(one_track_options[3]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(y_pattern_options.at(one_track_options[1]).y_tracks.at(one_track_options[3]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(y_pattern_options.at(one_track_options[1]).y_tracks.at(one_track_options[3]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.assumed_start_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
          temp_track.assumed_end_group = temp_track.assumed_start_group + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups() - 1;
          
          temp_track.track_line = one_track_line;
          
          temp_track.sort_cluster_lines_by_plate_index();
          
          temp_pattern.combined_tracks.push_back(temp_track);
          
          pattern_options.push_back(temp_pattern);
        }
        
        // cout << ">>>> " << option_best_distance[0] << endl;
      }
      else if (x_pattern_options.at(i).x_tracks.size() == 2 || y_pattern_options.at(j).y_tracks.size() == 2)
      {
        if (x_pattern_options.at(i).x_tracks.size() == 2)
        {
          for (int n = 0; n < y_pattern_options.at(j).y_tracks.size(); n++)
          {
            min_distance_on_diagonal_plates[0] = 0.;
            min_distance_on_diagonal_plates[1] = 0.;
            
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(0), y_pattern_options.at(j).y_tracks.at(n), temp_line);
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[k] = kk;
                }
              }
              
              min_distance_on_diagonal_plates[0] += best_dummy_distance;
            }
            
            if (min_distance_on_diagonal_plates[0] < option_best_distance[0])
            {
              option_best_distance[0] = min_distance_on_diagonal_plates[0];
              
              two_track_options[0][0] = i;
              two_track_options[0][1] = j;
              two_track_options[0][2] = 0;
              two_track_options[0][3] = n;
              two_track_options[0][4] = temp_closest_d_strip[0];
              two_track_options[0][5] = temp_closest_d_strip[1];
              two_track_line[0] = temp_line;
            }
            
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(1), y_pattern_options.at(j).y_tracks.at(n), temp_line);
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[k] = kk;
                }
              }
              
              min_distance_on_diagonal_plates[1] += best_dummy_distance;
            }
            
            if (min_distance_on_diagonal_plates[1] < option_best_distance[1])
            {
              option_best_distance[1] = min_distance_on_diagonal_plates[1];
              
              two_track_options[1][0] = i;
              two_track_options[1][1] = j;
              two_track_options[1][2] = 1;
              two_track_options[1][3] = n;
              two_track_options[1][4] = temp_closest_d_strip[0];
              two_track_options[1][5] = temp_closest_d_strip[1];
              two_track_line[1] = temp_line;
            }
          }
        }
        else if (y_pattern_options.at(j).y_tracks.size() == 2)
        {
          for (int m = 0; m < x_pattern_options.at(i).x_tracks.size(); m++)
          {
            min_distance_on_diagonal_plates[0] = 0.;
            min_distance_on_diagonal_plates[1] = 0.;
            
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(m), y_pattern_options.at(j).y_tracks.at(0), temp_line);
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[k] = kk;
                }
              }
              
              min_distance_on_diagonal_plates[0] += best_dummy_distance;
            }
            
            if (min_distance_on_diagonal_plates[0] < option_best_distance[0])
            {
              option_best_distance[0] = min_distance_on_diagonal_plates[0];
              
              two_track_options[0][0] = i;
              two_track_options[0][1] = j;
              two_track_options[0][2] = m;
              two_track_options[0][3] = 0;
              two_track_options[0][4] = temp_closest_d_strip[0];
              two_track_options[0][5] = temp_closest_d_strip[1];
              two_track_line[0] = temp_line;
            }
            
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(m), y_pattern_options.at(j).y_tracks.at(1), temp_line);
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[k] = kk;
                }
              }
              
              min_distance_on_diagonal_plates[1] += best_dummy_distance;
            }
            
            if (min_distance_on_diagonal_plates[1] < option_best_distance[1])
            {
              option_best_distance[1] = min_distance_on_diagonal_plates[1];
              
              two_track_options[1][0] = i;
              two_track_options[1][1] = j;
              two_track_options[1][2] = m;
              two_track_options[1][3] = 1;
              two_track_options[1][4] = temp_closest_d_strip[0];
              two_track_options[1][5] = temp_closest_d_strip[1];
              two_track_line[1] = temp_line;
            }
          }
        }
        
        if (option_best_distance[0] + option_best_distance[1] < 0.3) // !! use cut later
        {
          temp_pattern.reset();
          
          temp_pattern.pattern_quality[0] = option_best_distance[0] + option_best_distance[1];
          
          temp_track.cluster_lines.clear();
          temp_track.plate_index_of_cluster_line.clear();
          
          for (int k = 0; k < x_pattern_options.at(two_track_options[0][0]).x_tracks.at(two_track_options[0][2]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(x_pattern_options.at(two_track_options[0][0]).x_tracks.at(two_track_options[0][2]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(x_pattern_options.at(two_track_options[0][0]).x_tracks.at(two_track_options[0][2]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[0]].at(two_track_options[0][4]));
          temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[0]);
          if (number_of_d_plates > 1)
          {
            temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[1]].at(two_track_options[0][5]));
            temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[1]);
          }
          
          for (int k = 0; k < y_pattern_options.at(two_track_options[0][1]).y_tracks.at(two_track_options[0][3]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(y_pattern_options.at(two_track_options[0][1]).y_tracks.at(two_track_options[0][3]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(y_pattern_options.at(two_track_options[0][1]).y_tracks.at(two_track_options[0][3]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.assumed_start_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
          temp_track.assumed_end_group = temp_track.assumed_start_group + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups() - 1;
          
          temp_track.track_line = two_track_line[0];
          
          temp_track.sort_cluster_lines_by_plate_index();
          
          temp_pattern.combined_tracks.push_back(temp_track);
          
          temp_track.cluster_lines.clear();
          temp_track.plate_index_of_cluster_line.clear();
          
          for (int k = 0; k < x_pattern_options.at(two_track_options[1][0]).x_tracks.at(two_track_options[1][2]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(x_pattern_options.at(two_track_options[1][0]).x_tracks.at(two_track_options[1][2]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(x_pattern_options.at(two_track_options[1][0]).x_tracks.at(two_track_options[1][2]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[0]].at(two_track_options[1][4]));
          temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[0]);
          if (number_of_d_plates > 1)
          {
            temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[1]].at(two_track_options[1][5]));
            temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[1]);
          }
          
          for (int k = 0; k < y_pattern_options.at(two_track_options[1][1]).y_tracks.at(two_track_options[1][3]).cluster_lines.size(); k++)
          {
            temp_track.cluster_lines.push_back(y_pattern_options.at(two_track_options[1][1]).y_tracks.at(two_track_options[1][3]).cluster_lines.at(k));
            temp_track.plate_index_of_cluster_line.push_back(y_pattern_options.at(two_track_options[1][1]).y_tracks.at(two_track_options[1][3]).plate_index_of_cluster_line.at(k));
          }
          
          temp_track.assumed_start_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
          temp_track.assumed_end_group = temp_track.assumed_start_group + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups() - 1;
          
          temp_track.track_line = two_track_line[1];
          
          temp_track.sort_cluster_lines_by_plate_index();
          
          temp_pattern.combined_tracks.push_back(temp_track);
          
          pattern_options.push_back(temp_pattern);
        }
      }
      else
      {
        for (int m = 0; m < x_pattern_options.at(i).x_tracks.size(); m++)
        {
          for (int n = 0; n < y_pattern_options.at(j).y_tracks.size(); n++)
          {
            total_diagonal_distance = 0.;
            
            combine_3d_line_from_projections(x_pattern_options.at(i).x_tracks.at(m), y_pattern_options.at(j).y_tracks.at(n), temp_line);
            
            best_dummy_distance = 100.;
            
            for (int kk = 0; kk < cluster_lines[d_plate_indices[0]].size(); kk++)
            {
              find_line_to_line_distance(cluster_lines[d_plate_indices[0]].at(kk), temp_line, dummy_distance);
              
              if (dummy_distance < best_dummy_distance)
              {
                best_dummy_distance = dummy_distance;
                temp_closest_d_strip[0] = kk;
              }
            }
            
            total_diagonal_distance += best_dummy_distance;
            // cout << best_dummy_distance << " ";
            
            if (number_of_d_plates > 1)
            {
              best_dummy_distance = 100.;
              
              for (int kk = 0; kk < cluster_lines[d_plate_indices[1]].size(); kk++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[1]].at(kk), temp_line, dummy_distance);
                
                if (dummy_distance < best_dummy_distance)
                {
                  best_dummy_distance = dummy_distance;
                  temp_closest_d_strip[1] = kk;
                }
              }
              
              total_diagonal_distance += best_dummy_distance;
              // cout << best_dummy_distance << " ";
            }
            
            multi_track_best_distance[m][n] = total_diagonal_distance;
            multi_track_closest_d_strip[m][n][0] = temp_closest_d_strip[0];
            multi_track_closest_d_strip[m][n][1] = temp_closest_d_strip[1];
            multi_track_line[m][n] = temp_line;
            multi_track_line[m][n] = temp_line;
          }
          
          // cout << endl;
        }
        
        the_multi_track_selections.clear();
        operation_success = select_best_multi_track_options(multi_track_best_distance, x_pattern_options.at(i).x_tracks.size(), y_pattern_options.at(j).y_tracks.size(), the_multi_track_selections);
        
        /*
         if (operation_success)
         {
         cout << "true" << endl;
         }
         else
         {
         cout << "false" << endl;
         }
         */
        
        if (operation_success)
        {
          for (int t = 0; t < the_multi_track_selections.size(); t++)
          {
            //cout << "////// " << the_multi_track_selections.at(t).selection_quality << endl;
            //cout << the_multi_track_selections.at(t).index_x_track_selections[0] << " " << the_multi_track_selections.at(t).index_y_track_selections[0] << endl;
            //cout << the_multi_track_selections.at(t).index_x_track_selections[1] << " " << the_multi_track_selections.at(t).index_y_track_selections[1] << endl;
            //cout << the_multi_track_selections.at(t).index_x_track_selections[2] << " " << the_multi_track_selections.at(t).index_y_track_selections[2] << endl;
          }
        }
        
        for (int k = 0; k < the_multi_track_selections.size(); k++)
        {
          // if (the_multi_track_selections.at(k).selection_quality < 0.3) // !! make cut later
          {
            temp_pattern.reset();
            
            temp_pattern.pattern_quality[0] = the_multi_track_selections.at(k).selection_quality;
            
            for (int kk = 0; kk < the_multi_track_selections.at(k).number_of_tracks; kk++)
            {
              temp_track.cluster_lines.clear();
              temp_track.plate_index_of_cluster_line.clear();
              
              for (int s = 0; s < x_pattern_options.at(i).x_tracks.at(the_multi_track_selections.at(k).index_x_track_selections[kk]).cluster_lines.size(); s++)
              {
                temp_track.cluster_lines.push_back(x_pattern_options.at(i).x_tracks.at(the_multi_track_selections.at(k).index_x_track_selections[kk]).cluster_lines.at(s));
                temp_track.plate_index_of_cluster_line.push_back(x_pattern_options.at(i).x_tracks.at(the_multi_track_selections.at(k).index_x_track_selections[kk]).plate_index_of_cluster_line.at(s));
              }
              
              temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[0]].at(multi_track_closest_d_strip[the_multi_track_selections.at(k).index_x_track_selections[kk]][the_multi_track_selections.at(k).index_y_track_selections[kk]][0]));
              temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[0]);
              if (number_of_d_plates > 1)
              {
                temp_track.cluster_lines.push_back(cluster_lines[d_plate_indices[1]].at(multi_track_closest_d_strip[the_multi_track_selections.at(k).index_x_track_selections[kk]][the_multi_track_selections.at(k).index_y_track_selections[kk]][1]));
                temp_track.plate_index_of_cluster_line.push_back(d_plate_indices[1]);
              }
              
              for (int s = 0; s < y_pattern_options.at(j).y_tracks.at(the_multi_track_selections.at(k).index_y_track_selections[kk]).cluster_lines.size(); s++)
              {
                temp_track.cluster_lines.push_back(y_pattern_options.at(j).y_tracks.at(the_multi_track_selections.at(k).index_y_track_selections[kk]).cluster_lines.at(s));
                temp_track.plate_index_of_cluster_line.push_back(y_pattern_options.at(j).y_tracks.at(the_multi_track_selections.at(k).index_y_track_selections[kk]).plate_index_of_cluster_line.at(s));
              }
              
              temp_track.assumed_start_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
              temp_track.assumed_end_group = temp_track.assumed_start_group + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups() - 1;
              
              temp_track.track_line = multi_track_line[the_multi_track_selections.at(k).index_x_track_selections[kk]][the_multi_track_selections.at(k).index_y_track_selections[kk]];
              
              temp_track.sort_cluster_lines_by_plate_index();
              temp_pattern.combined_tracks.push_back(temp_track);
              
            }
            
            pattern_options.push_back(temp_pattern);
          }
        }
        
        
        // cout << "<<<<<<>>>>>>" << endl;
      }
      
      /*
       for (int k = 0; k < number_of_d_plates; k++)
       {
       best_dummy_distance = 100.;
       
       for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
       {
       find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line, dummy_distance);
       
       if (dummy_distance < best_dummy_distance)
       {
       best_dummy_distance = dummy_distance;
       }
       }
       
       min_distance_on_diagonal_plates[1] += best_dummy_distance;
       }
       */
      
      
    }
  }
  
  sort(pattern_options.begin(), pattern_options.end(), compare_pattern_quality_0);
  
  cout << "->->->-> " << pattern_options.size() << " ";
  for (int i = 0; i < pattern_options.size(); i++)
  {
    cout << pattern_options.at(i).pattern_quality[0] << " ";
  }
  cout << endl;
  
  if (pattern_options.size() > 0)
  {
    operation_success = true;
  }
  else
  {
    operation_success = false;
  }
  
  return operation_success;
}

bool track_finding_manager::select_best_multi_track_options(double multi_track_best_distance[][200], int number_of_x_tracks, int number_of_y_tracks, vector<multi_track_selections>& the_multi_track_selections)
{
  bool operation_success;
  
  double total_distance;
  
  multi_track_selections dummy_selection;
  
  if (number_of_x_tracks == 3 && number_of_y_tracks == 3)
  {
    total_distance = multi_track_best_distance[0][0] + multi_track_best_distance[1][1] + multi_track_best_distance[2][2];
// cout << total_distance << endl;
    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 0;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 1;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 2;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][1] + multi_track_best_distance[1][2] + multi_track_best_distance[2][0];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 1;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 2;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 0;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][2] + multi_track_best_distance[1][0] + multi_track_best_distance[2][1];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 2;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 0;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 1;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][2] + multi_track_best_distance[1][1] + multi_track_best_distance[2][0];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 2;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 1;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 0;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][1] + multi_track_best_distance[1][0] + multi_track_best_distance[2][2];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 1;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 0;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 2;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][0] + multi_track_best_distance[1][2] + multi_track_best_distance[2][1];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 0;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 2;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 1;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    sort(the_multi_track_selections.begin(), the_multi_track_selections.begin(), compare_multi_track_selections);
    
    operation_success = true;
  }
  else
  {
    operation_success = false;
  }
  
// cout << "++++++ " << the_multi_track_selections.size() << endl;
  
  return operation_success;
}

bool track_finding_manager::find_midstream_tracks_no_field(vector<pattern>& pattern_options, int event_class, track upstream_track)
{
  // It makes sense here to start all tracks at the vertex, so the point is common (as it really is),
  // which reduces the number of parameters. Therefore there are 3 (vertex position) + 3 x number_of_tracks
  // parameters in the fit. The vertex needs to have more weight, but doesn't it already?
  //
  
  
  
  
  return true;
}


bool track_finding_manager::propagate_tracks_downstream(vector<pattern>& pattern_options)
{
  // this function will eliminate the most divergent tracks that for sure will not reach downstream;
  // the rule here is if in doubt, leave it, it will not have clusters belonging to it anyway if
  // it doesn't actually propagate
  
  bool operation_success;
  
  plane_3d magnet_front_face;
  XYZVector intersection_point;
  
  double distance_to_axis; // distance between hit point on magnet front face plane and beam axis
  
  magnet_front_face.a = 0.;
  magnet_front_face.b = 0.;
  magnet_front_face.c = 1.;
  magnet_front_face.d = -8.; // magnet's front face is at 8 cm; do with file input later

  
  for (int i = 0; i < pattern_options.size(); i++)
  {
    cout << "tracks:" << pattern_options.at(i).combined_tracks.size() << endl;
    
    for (int j = 0; j < pattern_options.at(i).combined_tracks.size(); j++)
    {
      // cout << pattern_options.at(i).combined_tracks.at(j).track_line.point << " " << pattern_options.at(i).combined_tracks.at(j).track_line.direction << endl;
      
      operation_success = find_line_plate_intersection_point(pattern_options.at(i).combined_tracks.at(j).track_line, magnet_front_face, intersection_point);
      
      if (operation_success)
      {
        distance_to_axis = sqrt(intersection_point.X() * intersection_point.X() + intersection_point.Y() * intersection_point.Y());
        cout << distance_to_axis << " ";
        
        if (distance_to_axis > 3.) // !! fix later
        {
          pattern_options.at(i).combined_tracks.at(j).projects_downstream = false;
        }
        else
        {
          pattern_options.at(i).combined_tracks.at(j).projects_downstream = true;
        }
      }
      else
      {
        cout << "This shouldn't happen! Check the track_finding_manager::propagate_tracks_downstream code." << endl;
        cout << "Quitting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    cout << endl;
  }
  
  return operation_success;
}

bool track_finding_manager::assign_clusters_downstream_no_field(vector<pattern>& pattern_options, vector<line_3d> cluster_lines[])
{
  bool operation_success;
  bool geometric_calculation_success;
  
  int number_of_downstream_plates;
  int number_of_upstream_plus_midstream_plates;
  int number_of_total_plates;
  
  double dummy_distance;
  XYZVector intersection_point;
  ssd_plaque* p_dummy_plaque;
  
  number_of_downstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_plates();
  number_of_upstream_plus_midstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates() + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates();;
  number_of_total_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  // for now work with pattern option 0
  for (int i = 0; i < pattern_options.at(0).combined_tracks.size(); i++)
  {
    for (int j = number_of_upstream_plus_midstream_plates; j < number_of_total_plates; j++)
    {
      for (int k = 0; k < cluster_lines[j].size(); k++)
      {
        geometric_calculation_success = find_line_to_line_distance(pattern_options.at(0).combined_tracks.at(i).track_line, cluster_lines[j].at(k), dummy_distance);
        
        if (geometric_calculation_success && dummy_distance < 0.02) // !! do with cut; this is 200 microns
        {
          geometric_calculation_success = false;
          
          for (int l = 0; l < the_global_data_dispatcher->get_partial_geometry()->get_plate(j)->number_of_member_plaques; l++)
          {
            p_dummy_plaque = the_global_data_dispatcher->get_partial_geometry()->get_plaque(the_global_data_dispatcher->get_partial_geometry()->get_plate(j)->member_plaques[l]);
            
            if (find_line_plate_intersection_point(pattern_options.at(0).combined_tracks.at(i).track_line, p_dummy_plaque->plaque_plane_equation, intersection_point))
            {
              if (is_intersection_point_inside_plaque(intersection_point, p_dummy_plaque))
              {
                geometric_calculation_success = true;
              }
            }
          }
          
          if (geometric_calculation_success)
          {
            pattern_options.at(0).combined_tracks.at(i).cluster_lines.push_back(cluster_lines[j].at(k));
            pattern_options.at(0).combined_tracks.at(i).plate_index_of_cluster_line.push_back(j);
            
            // remove or don't remove cluster?
          }
        }
      }
    }
    
    pattern_options.at(0).combined_tracks.at(i).sort_cluster_lines_by_plate_index();
  }
  
  operation_success = true;
  
  return operation_success;
}

bool track_finding_manager::do_global_fit_no_field(vector<pattern>& pattern_options, track upstream_track)
{
  bool operation_success = true;
  
  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  double x_vertex, y_vertex, z_vertex;
  double dx_dz, dy_dz;
  
  // parameters:
  // 0: vertex.Z()
  // 1, 2 and later groups of two: direction of each track (px/pz and py/pz with dz = 1)
  // the number of parameters depends on the number of tracks
  
  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();

  initial_values_of_parameters.push_back(pattern_options[0].common_guessed_vertex.Z());

  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    initial_values_of_parameters.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction.X() / pattern_options[0].combined_tracks.at(i).track_line.direction.Z());
    initial_values_of_parameters.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction.Y() / pattern_options[0].combined_tracks.at(i).track_line.direction.Z());
  }
  
  errors_in_initial_values_of_parameters.push_back(0.2); // use input later
  
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    errors_in_initial_values_of_parameters.push_back(0.02); // ditto
    errors_in_initial_values_of_parameters.push_back(0.02);
  }
  
cout << "cccc" << endl;
  
  // the zeroth option should be the best; this must have been
  // sorted out by now, but there is still the possibility to use
  // any option - will be figured out with experimentation
  //
  my_global_fit_no_field_fcn the_fcn(pattern_options[0], upstream_track);
  
  VariableMetricMinimizer the_minimizer;
  
  FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
  cout << min.Fval() << endl;
  
  z_vertex = min.UserState().Value(0);
  
  x_vertex = upstream_track.track_line.point.X() + (z_vertex - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.X() / upstream_track.track_line.direction.Z();
  y_vertex = upstream_track.track_line.point.Y() + (z_vertex - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.Y() / upstream_track.track_line.direction.Z();
  
  pattern_options[0].vertex.SetCoordinates(x_vertex, y_vertex, z_vertex);
  
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    dx_dz = min.UserState().Value(i * 2 + 1);
    dy_dz = min.UserState().Value(i * 2 + 2);

    pattern_options[0].combined_tracks.at(i).track_line.direction.SetCoordinates(dx_dz, dy_dz, 1.);
    pattern_options[0].combined_tracks.at(i).track_line.direction = pattern_options[0].combined_tracks.at(i).track_line.direction.Unit();
  }
  
  pattern_options[0].pattern_quality[1] = min.Fval();
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->vertex = pattern_options[0].vertex;
  the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.clear();
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction);
  }
  the_global_data_dispatcher->get_event_reco_output_data_structure()->global_fit_quality = pattern_options[0].pattern_quality[1];
  
  for (int i = 0; i < pattern_options.size(); i++)
  {
    for (int j = 0; j < pattern_options.at(i).combined_tracks.size(); j++)
    {
      pattern_options.at(i).combined_tracks.at(j).find_intersection_points_rough();
    }
  }
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.clear();
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.push_back(pattern_options[0].combined_tracks.at(i));
  }
  
  return operation_success;
}
