module GenomeReader
  module CoordinateSystem
    # *.bed-format, UCSC internal format
    class ZeroBasedExclusive
      def seek_position(from)
        from < 0 ? 0 : from
      end

      def length(from, to)
        from = from < 0 ? 0 : from
        to - from
      end
    end

    # Ensembl, UCSC external(displaying) format, MAF
    class OneBasedInclusive
      def seek_position(from)
        (from < 1 ? 1 : from) - 1
      end

      def length(from, to)
        from = from < 1 ? 1 : from
        to - from + 1
      end
    end

    ZERO_BASED_EXCLUSIVE = ZeroBasedExclusive.new
    ONE_BASED_INCLUSIVE = OneBasedInclusive.new
  end
end
