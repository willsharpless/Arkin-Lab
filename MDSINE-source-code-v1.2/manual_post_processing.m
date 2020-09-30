% Manual Code for the Post Processing for Parameter.txt
config_filename = %'parameters_13.cfg'
config = formatconfig(config_filename);
disp('Performing post processing...')
post_processing(config)
disp(['Output written to ' config.general.output_dir])