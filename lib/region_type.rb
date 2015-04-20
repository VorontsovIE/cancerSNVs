require 'set'

class RegionType
  FEATURES = [:promoter, :intronic, :kataegis, :cage_peak, :in_repeat, :exon_coding]
  CALCULATED_FEATURES = [:regulatory]
  FEATURE_INQUIRIES = (FEATURES + CALCULATED_FEATURES).map{|feature| "#{feature}?".to_sym }

  attr_reader :features

  def initialize(features = Set.new)
    @features = features
  end

  def self.from_string(str)
    new( str.split(',').map(&:downcase).map(&:to_sym).to_set )
  end

  FEATURES.each do |feature|
    define_method "#{feature}?" do
      @features.include?(feature)
    end
  end

  def regulatory?
    intronic? || promoter?
  end


  def to_s
    @features.map(&:to_s).map(&:capitalize).join(',')
  end

  def <<(feature)
    @features << feature
  end

  def discard_types!
    @features = []
  end
end
