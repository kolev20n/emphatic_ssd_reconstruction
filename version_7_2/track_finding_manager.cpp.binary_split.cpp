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
  
  vector<track> midstream_partial_tracks[max_midstream_options];
  midstream_options the_midstream_options;
  int best_option_index;
  double best_option;
  
  double dummy_double;
  
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
      
      the_midstream_options.number_of_options = 0;
      
      if (operation_success)
      {
        operation_success = assign_clusters_midstream(cluster_lines, midstream_partial_tracks, the_midstream_options, event_class, upstream_track);
      }
      
      if (operation_success)
      {
        // for the moment work with the best option only; change this later maybe
        best_option_index = -1;
        best_option = -1.;
      
        for (int i = 0; i < the_midstream_options.number_of_options; i++)
        {
          if (the_midstream_options.quality_of_options[i] > best_option)
          {
            best_option_index = i;
            best_option = the_midstream_options.quality_of_options[i];
          }
        }
        
        operation_success = find_midstream_tracks_no_field(midstream_partial_tracks[best_option_index]);
      }
      
      // propagate tracks to downstream
      
      // assign clusters downstream
      
      // fit midstream plus downstream tracks (global fit?)
      
      // find vertex
      
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

bool track_finding_manager::assign_clusters_midstream(vector<line_3d> cluster_lines[], vector<track> midstream_partial_tracks[], midstream_options& the_midstream_options, int event_class, track upstream_track)
{
  // assuming the magnetic field is in the Y direction
  // the only four allowed directions here are X, Y and the two diagonal ones D1 and D2
  
  // need to keep track of unused clusters, so that further processing is possible (partial tracks etc.)
  
  // everywhere cluster_lines are considered in order; if not, it is a problem, because they are assigned in order to different tracks, so check !!
  
  bool operation_success = true;

  int number_of_x_plates = 0, number_of_y_plates = 0, number_of_d1_plates = 0, number_of_d2_plates = 0; // these are intended to measure X or Y respectively
  ssd_plate x_plates[max_number_of_plates_in_tracking_region], y_plates[max_number_of_plates_in_tracking_region], d1_plates[max_number_of_plates_in_tracking_region], d2_plates[max_number_of_plates_in_tracking_region];
  int x_plate_indices[max_number_of_plates_in_tracking_region],  y_plate_indices[max_number_of_plates_in_tracking_region], d1_plate_indices[max_number_of_plates_in_tracking_region], d2_plate_indices[max_number_of_plates_in_tracking_region];
  
  string plate_type;
  
  int u_plate, d_plate; // up- and downstream of the two plates
  int n_u_clusters, n_d_clusters; // number of u- and d-plate clusters
  
  vector<vector<track>> x_tracks, y_tracks;
  vector<track> temp_vector_track;
  track temp_track;
  
  bool success_y_scale = false, success_x_scale = false;
  
  double midstream_X_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_X_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  
  double scale_for_2;

  // move to cuts later; this is to account for misalignment without correction for the other direction
  double cut_min_tolerance = 0., cut_max_tolerance = 0.;

  double u_scale, d_scale;
  vector<double> projection_coordinates;
  
  bool decreasing_total_distance;
  XYZVector first_point, second_point;
  vector<double> first_distances, second_distances;
  double old_distance, new_distance;
  double start_z, mid_z, stop_z;
  double total_first_distance, total_second_distance;
  vector<double> first_projection_coordinates;
  vector<double> second_projection_coordinates;
  double dummy_double;
  double best_z;
  
//int scales_counter;
  
  for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates() + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates(); i++)
  {
    plate_type = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure;
    
    if (plate_type == "X")
    {
      x_plates[number_of_x_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      x_plate_indices[number_of_x_plates] = i;
      number_of_x_plates++;
    }
    else if (plate_type == "Y")
    {
      y_plates[number_of_y_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      y_plate_indices[number_of_y_plates] = i;
      number_of_y_plates++;
    }
    else if (plate_type == "D1")
    {
      d1_plates[number_of_d1_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      d1_plate_indices[number_of_d1_plates] = i;
      number_of_d1_plates++;
    }
    else if (plate_type == "D2")
    {
      d2_plates[number_of_d2_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      d2_plate_indices[number_of_d2_plates] = i;
      number_of_d2_plates++;
    }
    else
    {
      cout << "Unknown midstream plate type. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_X_scales(midstream_X_scales_max, midstream_X_scales_min);
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_Y_scales(midstream_Y_scales_max, midstream_Y_scales_min);
  
  // do for both directions; Y goes first because it is not affected by the magnetic field
  // and is therefore more important
  if (number_of_y_plates == 2)
  {
    u_plate = y_plate_indices[0];
    d_plate = y_plate_indices[1];
    
    n_u_clusters = cluster_lines[y_plate_indices[0]].size();
    n_d_clusters = cluster_lines[y_plate_indices[1]].size();
    
    if (n_u_clusters == 1 && n_d_clusters == 1)
    {
      y_tracks.clear();
      
      success_y_scale = true;
      
      y_tracks.clear();
      temp_vector_track.clear();
      
      temp_track.cluster_lines.clear();
      temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
      temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
      temp_vector_track.push_back(temp_track);
      y_tracks.push_back(temp_vector_track);
    }
    else if (n_u_clusters > 1 && n_d_clusters == 1)
    {
      y_tracks.clear();
      
      success_y_scale = true;
      
      for (int i = 0; i < n_u_clusters; i++)
      {
        temp_vector_track.clear();
        temp_track.cluster_lines.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_vector_track.push_back(temp_track);
        y_tracks.push_back(temp_vector_track);
      }
    }
    else if (n_u_clusters == 1 && n_d_clusters > 1)
    {
      y_tracks.clear();
      
      success_y_scale = true;
      
      for (int i = 0; i < n_d_clusters; i++)
      {
        temp_vector_track.clear();
        temp_track.cluster_lines.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
        temp_vector_track.push_back(temp_track);
        y_tracks.push_back(temp_vector_track);
      }
    }
    else if (n_u_clusters == 2 && n_d_clusters == 2)
    {
      //cout << "Y22" << endl;
      
      scale_for_2 = fabs((cluster_lines[d_plate].at(1).point.Y() - cluster_lines[d_plate].at(0).point.Y()) / (cluster_lines[u_plate].at(1).point.Y() - cluster_lines[u_plate].at(0).point.Y()));
      
      //cout << scale_for_2 << endl;
      
      // if more than 2 Y plates the scales are more than 1 number, so need to modify code
      if (scale_for_2 > midstream_Y_scales_min[0] - cut_min_tolerance && scale_for_2 < midstream_Y_scales_max[0] + cut_max_tolerance)
      {
        success_y_scale = true;
        
        y_tracks.clear();
        temp_vector_track.clear();

        temp_track.cluster_lines.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
        temp_vector_track.push_back(temp_track);
        
        temp_track.cluster_lines.clear();
        temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
        temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
        temp_vector_track.push_back(temp_track);
        
        y_tracks.push_back(temp_vector_track);
      }
    }
    else if (n_u_clusters > 2 && n_d_clusters == 2)
    {
      d_scale = fabs(cluster_lines[d_plate].at(1).point.Y() - cluster_lines[d_plate].at(0).point.Y());
      
      y_tracks.clear();
      
      for (int i = 0; i < n_u_clusters - 1; i++)
      {
        for (int j = n_u_clusters - 1; j > i; j--)
        {
          u_scale = fabs(cluster_lines[u_plate].at(j).point.Y() - cluster_lines[u_plate].at(i).point.Y());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_Y_scales_min[0] && d_scale / u_scale < midstream_Y_scales_max[0])
          {
            temp_vector_track.clear();
            
            temp_track.cluster_lines.clear();
            
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(i));
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(0));
            temp_vector_track.push_back(temp_track);
            
            temp_track.cluster_lines.clear();
            
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(j));
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(1));
            temp_vector_track.push_back(temp_track);
            
            y_tracks.push_back(temp_vector_track);
          }
        }
      }
      
      if (y_tracks.size() == 0)
      {
        success_y_scale = false;
      }
      else
      {
        success_y_scale = true;
      }
    }
    else if (n_u_clusters == 2 && n_d_clusters > 2)
    {
      u_scale = fabs(cluster_lines[u_plate].at(1).point.Y() - cluster_lines[u_plate].at(0).point.Y());
      
      y_tracks.clear();
      
      for (int i = 0; i < n_d_clusters - 1; i++)
      {
        for (int j = n_d_clusters - 1; j > i; j--)
        {
          d_scale = fabs(cluster_lines[d_plate].at(j).point.Y() - cluster_lines[d_plate].at(i).point.Y());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_Y_scales_min[0] && d_scale / u_scale < midstream_Y_scales_max[0])
          {
            temp_vector_track.clear();
            
            temp_track.cluster_lines.clear();
            
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(0));
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(i));
            temp_vector_track.push_back(temp_track);
            
            temp_track.cluster_lines.clear();
            
            temp_track.cluster_lines.push_back(cluster_lines[u_plate].at(1));
            temp_track.cluster_lines.push_back(cluster_lines[d_plate].at(j));
            temp_vector_track.push_back(temp_track);
            
            y_tracks.push_back(temp_vector_track);
          }
        }
      }
      
      if (y_tracks.size() == 0)
      {
        success_y_scale = false;
      }
      else
      {
        success_y_scale = true;
      }
    }
    else if (n_u_clusters > 2 && n_d_clusters > 2)
    {
      decreasing_total_distance = true;
     
      old_distance = 1000.;
      new_distance = 1000.;
      
      start_z = the_global_data_dispatcher->get_partial_geometry()->get_target_min_z();
      stop_z = the_global_data_dispatcher->get_partial_geometry()->get_target_max_z();
      mid_z = (stop_z + start_z) / 2.;
    
      while (decreasing_total_distance)
      {
        first_point = find_point_on_line_by_z(upstream_track.track_line, (start_z + mid_z) / 2.);
        second_point = find_point_on_line_by_z(upstream_track.track_line, (mid_z + stop_z) / 2.);

        first_projection_coordinates.clear();
        second_projection_coordinates.clear();

        for (int i = 0; i < n_u_clusters; i++)
        {
          dummy_double = fabs(first_point.Y() + (y_plates[0].position.Y() - first_point.Y()) * (y_plates[1].position.Z() - cluster_lines[u_plate].at(i).point.Z()) / (y_plates[0].position.Y() - first_point.Z()));
          
          if (dummy_double < 1.94) // do with plate size and cut (ssd_plate does not have a size)
          {
            first_projection_coordinates.push_back(dummy_double);
          }
          
          dummy_double = fabs(second_point.Y() + (y_plates[0].position.Y() - second_point.Y()) * (y_plates[1].position.Z() - cluster_lines[u_plate].at(i).point.Z()) / (y_plates[0].position.Y() - second_point.Z()));
          
          if (dummy_double < 1.94) // do with plate size and cut (ssd_plate does not have a size)
          {
            second_projection_coordinates.push_back(dummy_double);
          }
        }
        
        first_distances.clear();
        second_distances.clear();

        for (int i = 0; i < first_projection_coordinates.size(); i++)
        {
          first_distances.push_back(100.);
        }
        
        for (int i = 0; i < second_projection_coordinates.size(); i++)
        {
          second_distances.push_back(100.);
        }
        
        for (int i = 0; i < first_projection_coordinates.size(); i++)
        {
          for (int j = 0; j < n_d_clusters; j++)
          {
            if (fabs(first_projection_coordinates[i] - cluster_lines[d_plate].at(j).point.Y()) < first_distances.at(i))
            {
              first_distances.at(i) = fabs(first_projection_coordinates[i] - cluster_lines[d_plate].at(j).point.Y());
            }
          }
        }
        
        for (int i = 0; i < second_projection_coordinates.size(); i++)
        {
          for (int j = 0; j < n_d_clusters; j++)
          {
            if (fabs(second_projection_coordinates[i] - cluster_lines[d_plate].at(j).point.Y()) < second_distances.at(i))
            {
              second_distances.at(i) = fabs(second_projection_coordinates[i] - cluster_lines[d_plate].at(j).point.Y());
            }
          }
        }
        
        total_first_distance = 0.;
        total_second_distance = 0.;
        
        for (int i = 0; i < first_projection_coordinates.size(); i++)
        {
          total_first_distance += first_distances.at(i);
        }
        
        for (int i = 0; i < second_projection_coordinates.size(); i++)
        {
          total_second_distance += second_distances.at(i);
        }
        
        if (total_first_distance > total_second_distance)
        {
          new_distance = total_second_distance;
          start_z = mid_z;
          mid_z = (stop_z + start_z) / 2.;
        }
        else
        {
          new_distance = total_first_distance;
          stop_z = mid_z;
          mid_z = (stop_z + start_z) / 2.;
        }
        
        if (new_distance < old_distance)
        {
          old_distance = new_distance;
          best_z = mid_z;
        }
        else
        {
          decreasing_total_distance = false;
        }
cout << "---->" << new_distance << " " << best_z << endl;
      }
/*
scales_counter = 0;
    
      d_scale = fabs(cluster_lines[d_plate].at(0).point.Y() - cluster_lines[d_plate].at(n_d_clusters - 1).point.Y());
      
      for (int i = 0; i < n_u_clusters - 1; i++)
      {
        for (int j = n_u_clusters - 1; j > i; j--)
        {
          u_scale = fabs(cluster_lines[u_plate].at(j).point.Y() - cluster_lines[u_plate].at(i).point.Y());
          
          // add tolerances for alignment
          if (d_scale / u_scale > midstream_Y_scales_min[0] && d_scale / u_scale < midstream_Y_scales_max[0])
          {
cout << d_scale / u_scale << " ";
            scales_counter++;
          }
        }
      }
cout << endl;
cout << scales_counter << endl;
*/

      
    }
    else
    {
      cout << "Not enough clusters on a plate with presumed at least 1 cluster. " << endl;
      success_y_scale = false;
    }
  }
  else // more than 2 x or y planes midstream can be done later; not used in 2021 run
  {
    cout << "Midstream Y plates different than 2. Can't do cluster assignment. Quitting..." <<endl;
    exit(EXIT_FAILURE);
  }
  
  // now do the X direction
  if (number_of_x_plates == 2)
  {
    u_plate = x_plate_indices[0];
    d_plate = x_plate_indices[1];
    
    n_u_clusters = cluster_lines[x_plate_indices[0]].size();
    n_d_clusters = cluster_lines[x_plate_indices[1]].size();
  }
  else // more than 2 x or y planes midstream can be done later; not used in 2021 run
  {
    cout << "Midstream X plates different than 2. Can't do cluster assignment. Quitting..." <<endl;
    exit(EXIT_FAILURE);
  }
  
  return operation_success;
}

bool track_finding_manager::find_midstream_tracks_no_field(vector<track>& midstream_partial_tracks)
{
  // assuming the magnetic field is in the Y direction
  // only four directions are allowed: X, Y, and the two diagonal ones D1 and D2
  
  // need to keep track of unused clusters, so that further processing is possible (partial tracks etc.)
  
  // everywhere cluster_lines are considered in order; if not, it is a problem, because they are assigned in order to different tracks, so check !!
  
  bool operation_success;

  int number_of_x_plates = 0, number_of_y_plates = 0, number_of_d1_plates = 0, number_of_d2_plates = 0; // these are intended to measure X or Y respectively
  ssd_plate x_plates[max_number_of_plates_in_tracking_region], y_plates[max_number_of_plates_in_tracking_region], d1_plates[max_number_of_plates_in_tracking_region], d2_plates[max_number_of_plates_in_tracking_region];
  int  x_plate_indices[max_number_of_plates_in_tracking_region],  y_plate_indices[max_number_of_plates_in_tracking_region], d1_plate_indices[max_number_of_plates_in_tracking_region], d2_plate_indices[max_number_of_plates_in_tracking_region];
  
  string plate_type;

  for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates() + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates(); i++)
  {
    plate_type = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure;
    
    if (plate_type == "X")
    {
      x_plates[number_of_x_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      x_plate_indices[number_of_x_plates] = i;
      number_of_x_plates++;
    }
    else if (plate_type == "Y")
    {
      y_plates[number_of_y_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      y_plate_indices[number_of_y_plates] = i;
      number_of_y_plates++;
    }
    else if (plate_type == "D1")
    {
      d1_plates[number_of_d1_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      d1_plate_indices[number_of_d1_plates] = i;
      number_of_d1_plates++;
    }
    else if (plate_type == "D2")
    {
      d2_plates[number_of_d2_plates] = *the_global_data_dispatcher->get_partial_geometry()->get_plate(i);
      d2_plate_indices[number_of_d2_plates] = i;
      number_of_d2_plates++;
    }
    else
    {
      cout << "Unknown midstream plate type. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  return operation_success;
}
