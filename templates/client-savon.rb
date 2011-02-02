#!/usr/bin/env ruby

# Example SOAP client for the Mutalyzer webservice in Ruby using the savon
# library.
#
# See {path}/webservices
#
# Usage:
#   ruby client-savon.rb 'NM_002001.2:c.1del'
#
# This code is in the public domain; it can be used for whatever purpose
# with absolutely no restrictions.

require 'rubygems'
require 'savon'

URL = '{path}/services/?wsdl'

HTTPI.log = false

Savon.configure do |config|
  config.log = false
  config.log_level = :error
end

if ARGV.length == 0
  puts 'Please provide a variant'
  exit 1
end

client = Savon::Client.new do
  wsdl.document = URL
end

response = client.request :check_syntax do
  soap.body = {{ :variant => ARGV[0] }}
end

result = response.to_hash[:check_syntax_response][:check_syntax_result]

if result[:valid]
  puts 'Valid!'
else
  puts 'Not valid!'
end

if result[:messages]

  # This seems to be a bug in Savon. Arrays of length 1 are flattened,
  # so we cannot iterate over them.
  if not result[:messages][:soap_message].instance_of? Array
    result[:messages][:soap_message] = [result[:messages][:soap_message]]
  end

  result[:messages][:soap_message].each do |m|
    puts "Message (#{{m[:errorcode]}}): #{{m[:message]}}"
  end

end
