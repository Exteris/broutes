module Broutes
  module Maths
    class << self
      def haversine_distance(p1, p2)
        return unless p1.has_location? && p2.has_location?
        
        dlat = p2.lat - p1.lat
        dlon = p2.lon - p1.lon
    
        dlon_rad = dlon * RAD_PER_DEG
        dlat_rad = dlat * RAD_PER_DEG
    
        lat1_rad = p1.lat * RAD_PER_DEG
        lon1_rad = p1.lon * RAD_PER_DEG
    
        lat2_rad = p2.lat * RAD_PER_DEG
        lon2_rad = p2.lon * RAD_PER_DEG
    
        a = (Math.sin(dlat_rad/2))**2 + Math.cos(lat1_rad) * Math.cos(lat2_rad) * (Math.sin(dlon_rad/2))**2
        c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a))
    
        EARTH_RADIUS * c
      end

      # Calculate the median of arr
      # If the length is even, take the mean of the two medians
      def median(arr)
        sorted = arr.sort
        len = sorted.length
        m1 = sorted[(len - 1) / 2] 
        m2 = sorted[len / 2]
        if m1.nil? or m2.nil?
          return nil
        else
          return (m1 + m2) / 2.0
        end
      end

      def MAD(arr, m)
        median(arr.compact.map{|a| (a - m).abs})
      end
    end
  end
end
