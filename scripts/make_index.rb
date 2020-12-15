require 'nokogiri'

courses = ["cs653/src/", "phys516/src/", "cs596/src/"]
@content_children_map = {}

def make_index_for_directory(directory)
  children = @content_children_map[directory]

  fh = File.open(File.join(directory,"index.html"), "w")
  doc =  Nokogiri::HTML::Document.parse <<-EOHTML
  <h1>Index of #{ directory } </h1>
  <table>
  <tr><th>Name</th><th>Last modified</th><th>Size</th><th>Description</th></tr><tr><th colspan="5"><hr></th></tr>
  </table>
  EOHTML
  p doc.to_html
  body = doc.at_css "body"
  body.add_previous_sibling "<head>
  <title>Index of #{ directory } </title></head>"
  
  children.each do |child|
    tr            = doc.at_css "table"
    tr.add_child "<tr>
    <td><a href = #{ File.join(child) }> #{ child } </a></td><td>&nbsp;</td><td align=\"right\">  - </td><td>&nbsp;</td></tr>"
  end
  doc.write_to(fh)
end

def make_index(contents)
  @content_children_map = {}
  contents.each do |content| 
    children          = Dir.entries(content).select do |entry|
      !(entry == '.' || entry == '..' ) and !(entry == 'index.html')       
    end
    @content_children_map[content] = children   
  end
  p @content_children_map
  @content_children_map.each_key do |content|
    make_index_for_directory(content)
  end
end

make_index(courses)


