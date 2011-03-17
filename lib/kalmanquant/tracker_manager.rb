require 'set'
require 'ms/support/binary_search'


module Ms
  class Spectrum
    Peak = Struct.new(:spectrum, :index) do 
      def mz
        self[0].mzs[self[1]]
      end
      def intensity
        self[0].intensities[self[1]]
      end
    end
  end
end


module Kalmanquant
  class TrackerManager
    MAX_TRACKER_MISSES = 2
    MIN_TRACKER_SIZE = 4
    Match = Struct.new(:tracker, :peak, :distance)
    attr_accessor :spectra

    def initialize(spectra=[])
      @spectra = spectra
      @trackers = []
    end

    def find_features(max_misses=MAX_TRACKER_MISSES, min_size=MIN_TRACKER_SIZE)
      is_big_enough = lambda {|tracker| tracker.size >= min_size }

      active_trackers = []
      tracker_heaven = []
      #avg_delta_time = @spectra.map(:time).each_cons(2) {|a,b| b-a }.inject(:+) / @spectra.size-1
      #prev_time = @spectra.first.time - avg_delta_time
      @spectra.each do |spectrum|
        #delta_time = spectrum.time - prev_time
        (matched_trackers, unmatched_trackers, mz_index_set) = partition!(spectrum, active_trackers)
        (still_alive, dead) = unmatched_trackers.partition {|tracker| tracker.num_missed <= max_misses }
        (big_enough, died_miserly) = dead.partition &is_big_enough
        tracker_heaven.push(*big_enough)

        # create a new tracker for each unclaimed mz index
        new_trackers = []
        (0...spectrum.size).each do |mz_i|
          unless mz_index_set.include?(mz_i)
            new_trackers << Kalmanquant::Tracker.new(peak, 'mz_var_placeholder', 'int_var_placeholder', 'max_int_sqrt_placeholder') )
          end
        end
        # do new predictions for all active trackers
        (matched_trackers + still_alive).each {|tracker| tracker.predict! }
        active_trackers = matched_trackers + still_alive + new_trackers
      end
      (big_enough, too_small) = active_trackers.partition &is_big_enough
      tracker_heaven.push(*big_enough)
    end

    # greedy algorithm based on minimizing match distance and then negative
    # tracker.size.  returns a set of matching trackers, a list of unmatched
    # trackers, and the set of taken indices. As a side-effect, updates the
    # trackers info upon a match or not
    def partition!(spectrum, trackers)
      # also sort by confidence size (prefer smaller confidence intervals)
      sorted_matches = possible_matches(spectrum, trackers).sort_by {|match| [match.distance, match.tracker.confidence_size, -match.tracker.size] }
      taken_trackers = Set.new
      taken_peaks = Set.new
      sorted_matches.each do |match|
        unless taken_trackers.include?(match.tracker) || taken_peaks.include?(match.peak)
          peak = match.peak
          tracker = match.tracker
          tracker.update!(peak, spectrum.time)  # victory!
          taken_trackers.add(tracker)
          taken_peaks.add(peak)
        end
      end
      unmatched = trackers - taken_trackers
      unmatched.each {|tracker| tracker.missed! }
      [taken_trackers, unmatched, taken_peaks]
    end

    # returns an array of matches
    def possible_matches(spectrum, trackers)
      mzs = spectrum.mzs
      intensities = spectrum.intensities
      # brute force with binary search algorithm
      possible_matches = []
      trackers.each do |tracker|
        lo_mz_i = Ms::Support::BinarySearch.search_lower_boundary(mzs) {|mz| mz <=> tracker.lo_mz }
        # we use the range to ensure we only need to search the last part of the array
        if lo_mz_i
          hi_mz_i = Ms::Support::BinarySearch.search_upper_boundary(mzs, lo_mz_i...mzs.size) {|mz| mz <=> tracker.hi_mz }
          hi_mz_i -= 1 if hi_mz_i == mzs.size
        end
        if lo_mz_i && hi_mz_i
          (lo_mz_i..hi_mz_i).each do |i|
            if tracker.in_bounds?(mzs[i], intensities[i])
              possible_matches << Match.new(tracker, Ms::Spectrum::Peak(spectrum, i), tracker.normalized_distance(mzs[i], intensities[i]))
            end
          end
        end
      end
      possible_matches
    end
  end

end






=begin

classdef TrackerManager7 < handle
    %% Makes sure that the current and next scan data points are all
    % accounted for by the trackers. 
    % 
    properties (SetAccess = public)
        CurrentScanNum
        CurrentScan = [];
        ActiveTrackers = [];
        PICTrackers = [];
        
    end
    methods
        %% Constructor
        function trMgr = TrackerManager7(InitScanNum)
            %note that the scan num specified, may not 
            %be scan one
            trMgr.CurrentScanNum = InitScanNum;
        end
        function delete(h)
        end
        function [ActiveTrIndex PredPtsIndex Error] = predictActiveTrackers(trMgr)
            %Each tracker can only claim one point; therefore,
            %the max points claimed will be the number of 
            %active trackers
            NumActiveTrackers = size(trMgr.ActiveTrackers,1);
            PredPtsIndex = zeros(NumActiveTrackers, 1);
            ActiveTrIndex = PredPtsIndex;
            Error = PredPtsIndex;
            for i = 1:NumActiveTrackers
                %predict confidence intervals for Active Trackers
                predictCentroid(trMgr.ActiveTrackers(i));
                [DataPointIdx LeastError]= claimCurrentCentroid(trMgr.ActiveTrackers(i),...
                                                                trMgr.CurrentScan);
                %add only if claimed by tracker                             
                if DataPointIdx > 0                                 
                    PredPtsIndex(i,1) = DataPointIdx;
                    ActiveTrIndex(i,1) = i;
                    Error(i,1) = LeastError;
                end
            end
        end
        function manageStatus(trMgr, ActiveTrIndex, PredPtsIndex)
          
            %Make sure missed scans are accounted for
            manageMissed(trMgr, ActiveTrIndex == 0);
            %Call innovation step on all trackers that successfully claimed
            manageTracked(trMgr, PredPtsIndex(PredPtsIndex ~= 0), ActiveTrIndex ~= 0);
        end
        function manageMissed(trMgr, MissedIdx)
            MissedTrackers = trMgr.ActiveTrackers(MissedIdx);
            for i = 1:size(MissedTrackers,1)
                incrementMissed(MissedTrackers(i));
                assignStatus(MissedTrackers(i));
            end
        end
        function manageTracked(trMgr, DataPtsIndex, SelectedIndex)
                SelectedTrackers = trMgr.ActiveTrackers(SelectedIndex);
                for i = 1:size(SelectedTrackers,1)
                  
                   resetMissed(SelectedTrackers(i));
                   incrementTrackerLength(SelectedTrackers(i));
                   innovateCentroid(...
                            SelectedTrackers(i), ...
                            trMgr.CurrentScanNum, ...
                            DataPtsIndex(i), ...
                            trMgr.CurrentScan( DataPtsIndex(i), :)...
                            );     
                    assignStatus(SelectedTrackers(i));     
                end
                %Mark all tracked points unavailable to new trackers
                trMgr.CurrentScan(DataPtsIndex, :) = Inf;
        end
        function cleanActiveTrackers(trMgr, TrackerStatus)
            
            %initialize some arbitrarily large vectors to 
            %avoid dynamic array in loop
            maxScanPoints = 10^5;
            DeleteTrackerIndex = zeros(maxScanPoints,1);
            PICTrackersIndex = zeros(maxScanPoints,1);
            
            %do this so reshuffling only happens once
            for i = 1:size(trMgr.ActiveTrackers,1)
                if trMgr.ActiveTrackers(i).Status == TrackerStatus.NOISE
                    DeleteTrackerIndex(i) = i;
                elseif (trMgr.ActiveTrackers(i).Status == TrackerStatus.PIC) 
                    PICTrackersIndex(i) = i;
                end
            end
           
           %remember to remove zero elements from preallocation
           
           %Store Pure Ion Chromatogram (PIC)   
           %Consider removing peaks below a certain threshold.
           trMgr.PICTrackers = [trMgr.PICTrackers; ...
               trMgr.ActiveTrackers(PICTrackersIndex(PICTrackersIndex ~= 0), 1)];
           
           %Remove bad trackers and PIC trackers from list
           RmIndex = union(DeleteTrackerIndex(DeleteTrackerIndex ~= 0), ...
                           PICTrackersIndex);
           trMgr.ActiveTrackers = trMgr.ActiveTrackers(...
                                 setdiff(1:length(trMgr.ActiveTrackers),...
                                         RmIndex) ...
                                                   );
        end
        function initializeNewTrackers(trMgr, ghost, dt)
            %@param: 
            %       trMgr-> TrackerManager trMgrect
            %       ghost -> the initializer trMgrect
            %returns: nothing 
            
            %determine # of trackers to add
            ToClaim = trMgr.CurrentScan(:,2) ~= Inf;
            DataPoints = trMgr.CurrentScan(ToClaim,:);
            DataPointLength = size(DataPoints, 1);
            DataIndex = find(ToClaim);
            
            %no new trackers to add
            if (DataPointLength == 0)
               return; 
            end
            
            %every Tracker has common initialization
            Init = [ghost.MZVar, ghost.IntVar, ghost.MaxIntSqrt, trMgr.CurrentScanNum, dt];
               
			NewTrackers(DataPointLength,1) = Tracker7();
			trMgr.ActiveTrackers = [trMgr.ActiveTrackers; NewTrackers];
			   
            %every tracker has unique elements of initialization
            %id, mzseed, intseed, data num
            NewEnd =  size(trMgr.ActiveTrackers,1);
            NewStart = NewEnd - DataPointLength + 1;
            j = 1; %sub index for DataPointIndex
            for i = NewStart:NewEnd
                copyCommonInit(trMgr.ActiveTrackers(i, 1), Init);
                %all the following properties are dependent upon index: j
                trMgr.ActiveTrackers(i).MZXhat(1) = DataPoints(j, 1);
                trMgr.ActiveTrackers(i).IntXhat(1) = DataPoints(j, 2);
                trMgr.ActiveTrackers(i).DataPointNumbers(1) = DataIndex(j);
                trMgr.ActiveTrackers(i).ID = [trMgr.CurrentScanNum; DataIndex(j)]; 
                j = j + 1;
            end
        end
        function nonzerosPIC(trMgr)
           %clean up all the zeros in the trMgr at the end of the entire
           %scan
            for i = 1:size(trMgr.PICTrackers, 1)
                %% Data Structures to be kept in final implementation
                trMgr.PICTrackers(i).DataPointNumbers = ...
                    nonzeros(trMgr.PICTrackers(i).DataPointNumbers);
                trMgr.PICTrackers(i).ScanNumbers = ...
                    nonzeros(trMgr.PICTrackers(i).ScanNumbers);
              
                %% DEBUG MODE, temporary structures for figures
               trMgr.PICTrackers(i).MZPred = ...
                nonzeros(trMgr.PICTrackers(i).MZPred);
               trMgr.PICTrackers(i).IntPred = ...
                nonzeros(trMgr.PICTrackers(i).IntPred);                
                
                trMgr.PICTrackers(i).MZInno = ...
                    nonzeros(trMgr.PICTrackers(i).MZInno);
                trMgr.PICTrackers(i).IntInno = ...
                    nonzeros(trMgr.PICTrackers(i).IntInno);
                               
                trMgr.PICTrackers(i).TempPMZ = ...
                    nonzeros(trMgr.PICTrackers(i).TempPMZ);
                trMgr.PICTrackers(i).TempPInt = ...
                    nonzeros(trMgr.PICTrackers(i).TempPInt);
                
                trMgr.PICTrackers(i).InnoPMZ = ...
                    nonzeros(trMgr.PICTrackers(i).InnoPMZ);
                trMgr.PICTrackers(i).InnoPInt = ...
                    nonzeros(trMgr.PICTrackers(i).InnoPInt);
                
                tmpLeft = ...
                   nonzeros(trMgr.PICTrackers(i).MZConf(:,1));
                tmpRight = ...
                   nonzeros(trMgr.PICTrackers(i).MZConf(:,2));
                trMgr.PICTrackers(i).MZConf = [tmpLeft tmpRight];
                
                tmpUp = ...
                   nonzeros(trMgr.PICTrackers(i).IntConf(:,1));
                tmpDown = ...
                    nonzeros(trMgr.PICTrackers(i).IntConf(:,2));
                trMgr.PICTrackers(i).IntConf = [tmpUp tmpDown];      
            end
        end
    end
    methods (Static)
        %% Note I have never seen this function be used
        %Resolve Collisions
        function Map = competeActiveTrackers(Map, Error) 
             if TrackerManager7.existsCollision(nonzeros(Map)) == 1
                 Display('Resolving');
                 
                OneToOneMap = unique(Map);
                ItrEnd = size(OneToOneMap,1);
                if (OneToOneMap(1) == 0)
                   ItrEnd = Itr -1;
                end
                
                for i = 1:ItrEnd
                   %Get all indices of each unique point claimed
                   CollideIndex = find(Map == i);
                   %Two Trackers Claimed
                   if size(CollideIndex,1) > 1

                      %determine best tracker by least Error error in claim
                      [LeastError SubIndex] = min(Error(CollideIndex,1));
                      RejectSubIndex = find(CollideIndex ~= CollideIndex(SubIndex));

                      %make the datapoint inacessible to unsuccessful trackers
                      Map(CollideIndex(RejectSubIndex),1) = 0;
                   end
                end

                if TrackerManager7.existsCollision(Map) == 1
                    Display('Error in Resolving Collision: every centroid not resolved');
                end
             end
        end
        function TF = existsCollision(Map)
            %Is this a one-to-one mapping or not?
            %Assume it is not
            TF = 0;
            if size(Map, 2) ~= size(unique(Map),2)
                TF = 1;
                return;
            end
        end
    end
end %classdef

=end
