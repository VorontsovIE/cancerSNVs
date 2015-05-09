require 'set'

class RegionType
  FEATURES = [:promoter, :intronic, :kataegis]
  CALCULATED_FEATURES = [:regulatory]
  FEATURE_INQUIRIES = (FEATURES + CALCULATED_FEATURES).map{|feature| "#{feature}?".to_sym }

  attr_reader :features

  def initialize(features = Set.new)
    @features = features
  end

  def self.by_feature_status(promoter: false, intronic: false, kataegis: false)
    result = self.new.tap{|r|
      r << :kataegis  if kataegis
      r << :promoter  if promoter
      r << :intronic  if intronic
    }
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

  def ==(other)
    features == other.features
  end

  def equals?(other)
    other.is_a?(RegionType) &&features == other.features
  end

  def hash
    features.hash
  end

  def to_s
    @features.map(&:to_s).map(&:capitalize).join(',')
  end

  def description
    FEATURES.map{|feature|
      "#{feature}:#{features.include?(feature)}"
    }.join(' ')
  end

  def <<(feature)
    @features << feature
  end

  def discard_types!
    @features.clear
  end

  def self.each_possible
    return enum_for(:each_possible)  unless block_given?
    [true, false].repeated_permutation(3) do |promoter, intronic, kataegis|
      yield self.by_feature_status(promoter: promoter, intronic: intronic, kataegis: kataegis)
    end
  end
end
