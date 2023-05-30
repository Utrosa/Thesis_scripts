function fq=midi2freq(midi)
fq=(2.^((midi-69)/12))*440;
