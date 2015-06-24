def load_motif_underfitting_rates(fitting_log_filename)
  return {}  unless fitting_log_filename && File.exist?(fitting_log_filename)
  File.readlines(fitting_log_filename).drop(4).slice_before{|line|
    line.match(/^\t[^\t]/) && ! line.match(/(\d+) underfitted/) && ! line.match(/found \d+ from \d+/)
  }.map{|lines|
    motif = lines[0].strip
    underfitted = lines[1].strip.match(/^(\d+) underfitted/)[1].to_i
    [motif, underfitted]
  }.to_h
end
