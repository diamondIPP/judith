#include <iostream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>


int Clustered2RootTTree2Judith(const char* input_file_name, const char* output_file_name, const char* plane = "Plane0", const char* mode = "RECREATE", bool fill_event = true, Long64_t max_events = 0) {
	// this value depends on the input file
	static const UInt_t arr_size = 100000;

	// pyBAR ROOT file
	TFile* file = new TFile(input_file_name);
	TTree* table = (TTree*) file->Get("Table");

	// Judith ROOT file
	TFile* j_file = new TFile(output_file_name, mode);
// 	fill_event=false;   //bojan
	fill_event=true;
    if (fill_event)
		TTree* event = new TTree("Event", "Event");
	TDirectory* dir = j_file->mkdir(plane);
	dir->cd();
	TTree* hits = new TTree("Hits", "Hits");

	// event number and counter
	Long64_t envent_counter = 0;
	Long64_t curr_event_number = -1;

	// pyBAR arrays
	Int_t pybar_n_entries;
	Long64_t pybar_event_number[arr_size];
	UInt_t pybar_trigger_time_stamp[arr_size]; // bojan
	UShort_t pybar_cluster_ID[arr_size];
    UShort_t pybar_seed_column[arr_size];
    UShort_t pybar_seed_row[arr_size];
    Float_t pybar_mean_column[arr_size];
    Float_t pybar_mean_row[arr_size];
    UShort_t pybar_cluster_hits[arr_size];
    Float_t pybar_charge[arr_size];

	// Judith arrays
	Int_t judith_n_hits = 0;
    Long64_t judith_event_number=0;
// 	UInt_t judith_trigger_time_stamp; // bojan
	UShort_t judith_cluster_ID[arr_size];
	UShort_t judith_cluster_size[arr_size];
	Int_t judith_hit_pix_x[arr_size];
	Int_t judith_hit_pix_y[arr_size];
	Float_t judith_hit_pos_x[arr_size];
	Float_t judith_hit_pos_y[arr_size];
	UShort_t judith_hit_value[arr_size];

// 	ULong64_t judith_time_stamp;
	UInt_t judith_time_stamp;        //bojan
	ULong64_t judith_frame_number;
	Int_t judith_trigger_offset;
	Int_t judith_trigger_info;
	Bool_t judith_invalid;

	// pyBAR branches
	table->SetBranchAddress("n_entries", &pybar_n_entries);
	table->SetBranchAddress("event_number", pybar_event_number);
    table->SetBranchAddress("trigger_time_stamp", pybar_trigger_time_stamp);
    table->SetBranchAddress("ID", pybar_cluster_ID);
    table->SetBranchAddress("n_hits", pybar_cluster_hits);
	table->SetBranchAddress("charge", pybar_charge);
	table->SetBranchAddress("seed_column", pybar_seed_column);
	table->SetBranchAddress("seed_row", pybar_seed_row);
	table->SetBranchAddress("mean_column", pybar_mean_column);
	table->SetBranchAddress("mean_row", pybar_mean_row);
// 	table->SetBranchAddress("trigger_time_stamp", pybar_trigger_time_stamp);
// 	table->SetBranchAddress("trigger_number", pybar_trigger_time_stamp);
// 	table->SetBranchAddress("relative_BCID", pybar_relative_bcid);
// 	table->SetBranchAddress("event_status", pybar_event_status);

	// Judith branches
	hits->Branch("Event_number", &judith_event_number, "Event_number/L");
// 	hits->Branch("TimeStamp", &judith_trigger_time_stamp, "TimeStamp/I");
	hits->Branch("NClusters", &judith_n_hits, "NClusters/I");
	hits->Branch("ID", judith_cluster_ID, "ClusID[NClusters]/s");
	hits->Branch("ClusterSize", judith_cluster_size, "ClusSize[NClusters]/s");
	hits->Branch("PixX", judith_hit_pix_x, "HitPixX[NClusters]/I");
	hits->Branch("PixY", judith_hit_pix_y, "HitPixY[NClusters]/I");
	hits->Branch("PosX", judith_hit_pos_x, "HitPosX[NClusters]/F");
	hits->Branch("PosY", judith_hit_pos_y, "HitPosY[NClusters]/F");
	hits->Branch("Value", judith_hit_value, "HitValue[NClusters]/s");
	
	if (fill_event) {
      event->Branch("TimeStamp", &judith_time_stamp, "TimeStamp/I");
//       event->Branch("FrameNumber", &judith_frame_number, "FrameNumber/l");
//       event->Branch("TriggerOffset", &judith_trigger_offset, "TriggerOffset/I");
//       event->Branch("TriggerInfo", &judith_trigger_info, "TriggerInfo/I");
//       event->Branch("Invalid", &judith_invalid, "Invalid/O");
	}
	const Long64_t max_chunks = table->GetEntries();

	for (Long64_t curr_chunk = 0; curr_chunk < max_chunks; curr_chunk++) {
		if (max_events > 0 && envent_counter >= max_events - 1) {
          break;
		}
		int chunk_bytes = (int) table->GetEntry(curr_chunk);
                
//                 cout << pybar_trigger_time_stamp[1] << endl;

		int chunk_entries = (int) pybar_n_entries;
		std::cout << "reading chunk " << curr_chunk << " with size "
				<< chunk_entries << std::endl;

		for (int curr_chunk_index = 0; curr_chunk_index < chunk_entries;
				curr_chunk_index++) {
			// new event
			if (pybar_event_number[curr_chunk_index] != curr_event_number) {
				// in case max_events is set and reached
				if (max_events > 0 && envent_counter >= max_events - 1) {
					std::cout << "reached max. events " << max_events << " at chunk "
							<< curr_chunk << " index " << curr_chunk_index
							<< std::endl;
					break;
				}
				// store_event is set to false during initialization to prevent writing empty event
				if (!(curr_chunk == 0 && curr_chunk_index == 0)) {
					hits->Fill();
					if (fill_event)
						event->Fill();
					envent_counter++;
				}
				curr_event_number = pybar_event_number[curr_chunk_index];
				judith_n_hits = 0;
                judith_event_number++;
				// fill event TTree with new event
				if (fill_event) {
					// fill event
// 					judith_time_stamp = (ULong64_t) pybar_trigger_time_stamp[curr_chunk_index];
//                                         judith_frame_number = (ULong64_t) pybar_event_number[curr_chunk_index];
// 					judith_trigger_offset = 0;
// 					judith_trigger_info = 0;
// 					// check for unknown words
// 					if ((pybar_event_status[curr_chunk_index] & 0b0000000000010000)
// 							== 0b0000000000010000) {
// 						judith_invalid = 1;
// 					} else {
// 						judith_invalid = 0;
// 					}
                  judith_time_stamp = pybar_trigger_time_stamp[0];
//                   cout << pybar_trigger_time_stamp[0] << " " << judith_time_stamp << endl;

				}
			// in case store_event is false loop over remaining hits of the current event
			}


			// empty event
			if (0) {
				;
			} 
			else if (judith_n_hits >= arr_size) {
				std::cout << "reached the array size limit at chunk "
						<< curr_chunk << " index " << curr_chunk_index
						<< "event" << curr_event_number
						<< std::endl;
				;
			} 
			else {
				// fill hits TTree
				// pyBAR: starting col / row from 1
				// Judith: starting col / row from 0
// 				judith_event_number = pybar_event_number;
                judith_cluster_ID[judith_n_hits] = pybar_cluster_ID[curr_chunk_index];
                judith_cluster_size[judith_n_hits] = pybar_cluster_hits[curr_chunk_index];
                judith_hit_pix_x[judith_n_hits] = pybar_seed_column[curr_chunk_index] - 1;
				judith_hit_pix_y[judith_n_hits] = pybar_seed_row[curr_chunk_index] - 1;
				judith_hit_value[judith_n_hits] = (UShort_t) pybar_charge[curr_chunk_index];
// 				judith_hit_timing[judith_n_hits] = (Int_t) pybar_relative_bcid[curr_chunk_index];
// 				judith_hit_timing[judith_n_hits] = 0;
// 				judith_hit_in_cluster[judith_n_hits] = -1;
				judith_hit_pos_x[judith_n_hits] = pybar_mean_column[curr_chunk_index];
				judith_hit_pos_y[judith_n_hits] = pybar_mean_row[curr_chunk_index];
// 				judith_hit_pos_z[judith_n_hits] = 0;
				judith_n_hits++;
//                 cout << pybar_cluster_ID << " " << pybar_cluster_hits << endl;
//                 cout << judith_cluster_ID[judith_n_hits] << " " << judith_cluster_size[judith_n_hits] << endl;
                
                
			}
			// reached end of chunk and leave a message
			if (curr_chunk_index == chunk_entries - 1) {
				std::cout << "reached the end of chunk at chunk "
						<< curr_chunk << " index " << curr_chunk_index
						<< std::endl;
			}
		}
		// reached max_chunks and leave a message
		if (curr_chunk == max_chunks - 1) {
			std::cout << "reached the end of file at chunk "
					<< curr_chunk << std::endl;
		}
	}

	hits->Fill();
	if (fill_event)
		event->Fill();
	envent_counter++;

	j_file->Write();
	j_file->Close();
    
    cout << judith_event_number << endl;
    cout << envent_counter << endl;
    gApplication->Terminate();
    return envent_counter;
}
