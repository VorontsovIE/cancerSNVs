$:.unshift File.absolute_path('lib', __dir__)
$:.unshift '/home/ilya/iogen/projects/cage_analysis/lib/intervals/' # File.absolute_path('lib', __dir__)
require 'genome_region'

# class IntervalTree
#   attr_reader :chromosomes
#   def initialize(interval_object_pairs)
#     @chromosomes = interval_object_pairs
#       .group_by{|interval, object| "#{interval.chromosome},#{interval.strand}" }
#       .map{|chr, interval_object_pairs_subset|
#         # p interval_object_pairs_subset
#         intervals = interval_object_pairs_subset
#         .each_with_index
#         .sort_by{|(interval, object), ind|
#           interval.pos_end
#         }
#         [chr,intervals]
#       }.to_h#.tap{|x| p x}
#   end

#   def find(chr, strand, pos)
#     (interval, obj), ind = @chromosomes["#{chr},#{strand}"].bsearch{|(interval, obj), ind| interval.pos_end >= pos } # ?
#     @chromosomes["#{chr},#{strand}"].drop(ind).take_while{|(interval, obj), ind| interval.pos_end}
#     interval.include_position?(pos) ? [interval, obj] : nil
#     # return []  unless ind
#     # # select tss-es in upstream and in downstream (because promoter region can be assymetric)
#     # tsses[ [ind - 1, 0].max, 2 ].map{|tss,ind| tss}
#   end
# end
#
# interval_tree = IntervalTree.new(interval_object_pairs.first(20))
# p interval_tree.find('CHR_HSCHR6_MHC_QBL_CTG1', :+, 29887830)
