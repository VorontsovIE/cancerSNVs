class Table
  # Absent header null-object
  class EmptyHeader
    attr_accessor :size
    
    def initialize(size)
      @size = size
    end
    
    def header_names
      Array.new(size){ [nil, ''] }
    end
    
    def to_s
      ''
    end

    def inspect
      'EmptyHeader'
    end
    
    def header_index(header_symbol)
      nil
    end
    
    def header_name(header_symbol)
      nil
    end
    
    def etc(*args)
      'empty'
    end
    
    def empty?
      true
    end
    
    def +(other)
      other.empty? ? EmptyHeader.new(size + other.size) : Header.new(header_names + other.header_names)
    end
    
    def ==(other)
      other.empty? && size == other.size
    end

    def each
      return enum_for(:each)  unless block_given?
      header_names.each{|(k,v)| yield v }
    end
    include Enumerable

    def select_columns(columns)
      EmptyHeader.new(columns.size)
    end
  end

  # Header (it can be used as a row or as a column header)
  class Header
    def self.create(headers)
      return headers  if headers.is_a?(Header) || headers.is_a?(EmptyHeader)
      return Header.new([[(headers.empty? ? nil : headers), headers],])  if headers.is_a?(String) || headers.is_a?(Symbol)
      if headers.is_a?(Enumerable)
        header_names = headers.map do |header|
          if header.is_a?(Array) && header.size == 2
            header
          elsif header.is_a?(String) || header.is_a?(Symbol)
            [header, header]
          else
            raise "Unknown header format: #{headers.inspect}"
          end
        end
        Header.new(header_names)
      else
        raise "Unknown header format: #{headers.inspect}"
      end
    end

    def self.common_header(set_of_headers)
      non_empty_headers = set_of_headers.select{|hdr| ! hdr.empty? }
      suitable = non_empty_headers.size < 2 || non_empty_headers.all?{|hdr| hdr == non_empty_headers.first }
      raise "Non-consistent header columns:\n #{set_of_headers.map(&:inspect).join("\n")}"  unless suitable
      non_empty_headers.first
    end

    attr_accessor :header_names # [[:col_1, 'column 1 name'], [:col_2, 'column 2 name'], ...]

    def initialize(header_names)
      @header_names = header_names
      symbols = @header_names.map{|k,v| k }.compact
      raise "Duplicated header symbols in: #{symbols.inspect}"  unless symbols.uniq.size == symbols.size
    end

    def header_index(header_symbol)
      @header_names.index{|k,v| k == header_symbol }
    end
    
    def header_name(header_symbol)
      @header_name[ header_index(header_symbol) ]
    end

    def to_s
      @header_names.map{|k,v| v }.join("\t")
    end

    def inspect
      "TableHeader<#{etc}>"
    end

    def size
      @header_names.size
    end

    def etc(first_n: 5, separator: '; ', etc_sign: '; ...')
      @header_names.first(first_n).map{|k,v| v }.join(separator) + (size > first_n ? etc_sign : '')
    end

    def empty?
      false
    end

    def +(other)
      Header.new(header_names + other.header_names)
    end

    def ==(other)
      header_names == other.header_names
    end

    def each
      return enum_for(:each)  unless block_given?
      header_names.each{|(k,v)| yield v }
    end
    include Enumerable

    def select_columns(columns)
      col_indices = columns.map{|column| header_index(column) }
      Header.new(col_indices.map{|col_index| @header_names[col_index]})
    end
  end
end
