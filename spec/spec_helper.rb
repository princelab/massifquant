require 'rubygems'
require 'spec/more'

$LOAD_PATH.unshift(File.dirname(__FILE__))
$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
require 'kalmanquant'

Bacon.summary_on_exit
