require 'bioinform'

raise 'Specify folder with collection'  unless from_folder = ARGV[0] # '/home/ilya/iogen/hocomoco/'
raise 'Specify output folder'  unless output_folder = ARGV[1] # './source_data/hocomoco_expanded_flanks'

Dir.mkdir(output_folder)  unless Dir.exist?(output_folder)
Dir.glob(File.join(from_folder, '*.pwm')) {|fn|
  pwm = Bioinform::MotifModel::PWM.from_file(fn)
  # pwm_expanded = pwm.left_augmented(pwm.length / 2).right_augmented(pwm.length - pwm.length / 2)
  pwm_expanded = pwm.left_augmented(11).right_augmented(11)
  File.write(File.join(output_folder, File.basename(fn)), pwm_expanded.to_s)
}
