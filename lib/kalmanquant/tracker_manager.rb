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
    KEEP_DEAD_RUNTS = true
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
      new_trackers = []
      tracker_heaven = []
      trackers_died_small = []

      @spectra.each do |spectrum|
        active_trackers.each {|tracker| tracker.predict!(spectrum.time) }

        (matched_trackers, unmatched_trackers, matched_peaks, unmatched_peaks) = compete(spectrum, active_trackers + new_trackers)
        unmatched_trackers.each {|tracker| tracker.missed! }

        # send dead trackers to trash or tracker heaven
        (still_alive, dead) = unmatched_trackers.partition {|tracker| tracker.num_missed <= max_misses }
        # we divide out the dead in case we don't want to keep little trackers
        # around out of memory concerns
        (big_enough, not_big_enough) = dead.partition &is_big_enough
        tracker_heaven.push(*big_enough) 
        trackers_died_small.push(*not_big_enough) if KEEP_DEAD_RUNTS

        # create a new tracker for each unclaimed mz index
        new_trackers = unmatched_peaks.map do |peak|
          Kalmanquant::Tracker.new(peak, spectrum.time, ghost_peak)
        end

        matched_trackers.zip(matched_peaks) do |tracker,peak|
          tracker.update!(peak, spectrum.time)
        end

        active_trackers = matched_trackers + still_alive
      end
      # 2nd coming ... judgement for everyone! hurray
      (big_enough, too_small) = (active_trackers+new_trackers).partition &is_big_enough
      tracker_heaven.push(*big_enough)
    end

    def compete(spectrum, trackers)
      sorted_matches = possible_matches(spectrum, trackers).sort_by {|match| [match.distance, match.tracker.confidence_size, -match.tracker.size] }
      set_of_matched_trackers = Set.new
      set_of_matched_peaks = Set.new
      matched_trackers, matched_peaks, unmatched_trackers, unmatched_peaks = [], [], [], []
      sorted_matches.each do |match| 
        tracker, peak = match.tracker, match.peak
        if set_of_matched_trackers.include?(tracker) || set_of_matched_peaks.include?(peak)
          unmatched_trackers << tracker && unmatched_peaks << peak
        else
          set_of_matched_trackers << tracker && set_of_matched_peaks << peak
          matched_trackers << tracker && matched_peaks << peak
        end
      end
      [matched_trackers, unmatched_trackers, matched_peaks, unmatched_peaks]
    end

    # returns an array of matches
    def possible_matches(spectrum, trackers)
      mzs = spectrum.mzs
      intensities = spectrum.intensities
      # brute force with binary search algorithm
      possible_matches = []
      trackers.each do |tracker|
        (lo_mz, hi_mz) = tracker.pred_mz_margins
        lo_mz_i = Ms::Support::BinarySearch.search_lower_boundary(mzs) {|mz| mz <=> lo_mz }
        # we use the range to ensure we only need to search the last part of the array
        if lo_mz_i
          hi_mz_i = Ms::Support::BinarySearch.search_upper_boundary(mzs, lo_mz_i...mzs.size) {|mz| mz <=> hi_mz }
          hi_mz_i -= 1 if hi_mz_i == mzs.size
        end
        if lo_mz_i && hi_mz_i
          (lo_mz_i..hi_mz_i).each do |i|
            # the in_bounds? could be asked by checking if the normalized
            # distance is within desired bounds ....!!!!
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

