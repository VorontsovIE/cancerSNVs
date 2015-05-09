require 'set'

class RegionType
  FEATURES = [:promoter, :intronic, :kataegis]
  CALCULATED_FEATURES = [:regulatory]
  FEATURE_INQUIRIES = (FEATURES + CALCULATED_FEATURES).map{|feature| "#{feature}?".to_sym }

  attr_reader :features

  def initialize(features = Set.new)
    @features = features
  end

  def self.from_string(str)
    # # Optimized version of:
    # new( str.split(',').lazy.map(&:downcase).map(&:to_sym).to_set )
    features = Set.new
    str.split(',').each do |feature|
      features << feature.downcase.to_sym
    end
    new(features)
  end

  FEATURES.each do |feature|
    define_method "#{feature}?" do
      @features.include?(feature)
    end
  end

  def regulatory?
    (intronic? || promoter?) && ! kataegis?
  end


  def to_s
    @features.map(&:to_s).map(&:capitalize).join(',')
  end

  def <<(feature)
    @features << feature
  end

  def discard_types!
    @features.clear
  end
end
