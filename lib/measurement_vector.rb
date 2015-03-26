class Vector
  attr_reader :values
  def initialize(values)
    @values = values
  end

  def [](ind)
    @values[ind]
  end

  def dim
    values.size
  end

  def +(other)
    raise 'Incompatible dimensionality'  unless dim == other.dim
    Vector.new( dim.times.map{|i| self[i] + other[i]} )
  end

  def -(other)
    raise 'Incompatible dimensionality'  unless dim == other.dim
    Vector.new( dim.times.map{|i| self[i] - other[i]} )
  end

  def *(other)
    if other.is_a?(Vector)
      raise 'Incompatible dimensionality'  unless dim == other.dim
      Vector.new( dim.times.map{|i| self[i] * other[i] } )
    elsif other.is_a?(Numeric)
      Vector.new( @values.map{|el| el * other } )
    else
      raise
    end
  end

  def /(other)
    if other.is_a?(Vector)
      raise 'Incompatible dimensionality'  unless dim == other.dim
      Vector.new( dim.times.map{|i| self[i] / other[i] } )
    elsif other.is_a?(Numeric)
      Vector.new( @values.map{|el| el / other } )
    else
      raise
    end
  end

  def **(k)
    Vector.new( @values.map{|el| el ** k } )
  end

  def coerce(other)
    if other.is_a?(Numeric)
      [self, other]
    else
      raise
    end
  end

  def each(&block)
    @values.each(&block)
  end

  def to_s
    @values.join("\t")
  end

  def inspect
    '<' + @values.join(',') + '>'
  end

  def round(rate = 0)
    Vector.new(@values.map{|el| el.round(rate) })
  end

  include Enumerable
end
