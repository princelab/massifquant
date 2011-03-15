
module Kalmanquant
  # features are detected in reverse order
  # returns lc peaks
  def self.find_features(array_of_spectra, reverse_order=true)
    array_of_spectra.reverse! if reverse_order



    array_of_spectra.reverse! if reverse_order
    lc_peaks
  end
end
