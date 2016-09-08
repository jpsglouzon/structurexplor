$(function() {    
    
    /////////Export forna structure in svg
    $("#export-rna_ss1-svg").on ("click", function (e) {

        datamonkey.save_image("svg", "#rna_ss1");
        ga('send', 'event', 'widget', 'Export structure 1 in SVG');


    });
    $("#export-rna_ss2-svg").on ("click", function (e) {

        datamonkey.save_image("svg", "#rna_ss2");
       ga('send', 'event', 'widget', 'Export structure 1 in SVG');

    });
    
    $("#export-phylo-svg").on ("click", function (e) {

        datamonkey.save_image("svg", "#dendSS2");
        ga('send', 'event', 'widget', 'Export hierarchy in SVG');

    });
    
    $("#export-scatplotsnm-svg").on ("click", function (e) {

      datamonkey.save_image("svg", "#scatplotsnm");
      ga('send', 'event', 'widget', 'Export snm scatterplot in SVG');


    });
    
    $("#export-varexplained-svg").on ("click", function (e) {

    datamonkey.save_image("svg", "#varExp");
          ga('send', 'event', 'widget', 'Export explained variability');


    });
    
    $('#fullScreen').click(function() {
      $('#step72').css({
          "zindex": "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999", 
          "width": "100%",
          "height": "100%", 
          "position": "relative", 
          "top": "0", 
          "left": "0", 
           "background-color": "white"

      });
   });
    
    
      var root = this;
      
      var datamonkey = function () {};
      
      if (typeof exports !== 'undefined') {
        if (typeof module !== 'undefined' && module.exports) {
          exports = module.exports = Datamonkey;
        }
        exports.datamonkey = datamonkey;
      } else {
        root.datamonkey = datamonkey;
      }
          
      datamonkey.save_image = function(type, container) {
      
        var prefix = {
          xmlns: "http://www.w3.org/2000/xmlns/",
          xlink: "http://www.w3.org/1999/xlink",
          svg: "http://www.w3.org/2000/svg"
        }
      
        function get_styles(doc) {
      
          function process_stylesheet(ss) {
            try {
              if (ss.cssRules) {
                for (var i = 0; i < ss.cssRules.length; i++) {
                  var rule = ss.cssRules[i];
                  if (rule.type === 3) {
                    // Import Rule
                    process_stylesheet(rule.styleSheet);
                  } else {
                    // hack for illustrator crashing on descendent selectors
                    if (rule.selectorText) {
                      if (rule.selectorText.indexOf(">") === -1) {
                        styles += "\n" + rule.cssText;
                      }
                    }
                  }
                }
              }
            } catch (e) {
              console.log('Could not process stylesheet : ' + ss);
            }
          }
      
         // var styles = "",
          var    styleSheets = doc.styleSheets;
      
          if (styleSheets) {
            for (var i = 0; i < styleSheets.length; i++) {
              process_stylesheet(styleSheets[i]);
            }
          }
      
          return styles;
      
        }
      
        var convert_svg_to_png = function(image_string) {
      
          var image = document.getElementById("hyphy-chart-image");
      
          image.onload = function() {
      
            var canvas = document.getElementById("hyphy-chart-canvas");
            canvas.width = image.width;
            canvas.height = image.height;
            var context = canvas.getContext("2d");
            context.fillStyle = "#FFFFFF";
            context.fillRect(0,0,image.width,image.height);
            context.drawImage(image, 0, 0);
            var img = canvas.toDataURL("image/png");
            var pom = document.createElement('a');
            pom.setAttribute('download', 'image.png');
            pom.href = canvas.toDataURL("image/png");     
            $("body").append(pom);
            pom.click();
            pom.remove();
      
          }
      
          image.src = image_string;
      
        }
      
        var svg = $(container).find("svg")[0];
        if (!svg) {
          svg = $(container)[0];
        }
      
        var styles = get_styles(window.document);
      
        svg.setAttribute("version", "1.1");
      
        var defsEl = document.createElement("defs");
        svg.insertBefore(defsEl, svg.firstChild); 
      
        var styleEl = document.createElement("style")
        defsEl.appendChild(styleEl);
        styleEl.setAttribute("type", "text/css");
      
        // removing attributes so they aren't doubled up
        svg.removeAttribute("xmlns");
        svg.removeAttribute("xlink");
      
        // These are needed for the svg
        if (!svg.hasAttributeNS(prefix.xmlns, "xmlns")) {
          svg.setAttributeNS(prefix.xmlns, "xmlns", prefix.svg);
        }
      
        if (!svg.hasAttributeNS(prefix.xmlns, "xmlns:xlink")) {
          svg.setAttributeNS(prefix.xmlns, "xmlns:xlink", prefix.xlink);
        }
      
        var source = (new XMLSerializer()).serializeToString(svg).replace('</style>', '<![CDATA[' + styles + ']]></style>').replace('<defs xmlns="http://www.w3.org/1999/xhtml">', '<defs>').replace('</defs><style/>', '</defs>'        );
        var rect = svg.getBoundingClientRect();
        var doctype = '<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';
        var to_download = [doctype + source]
        var image_string = 'data:image/svg+xml;base66,' + encodeURIComponent(to_download);
      
        if(type == "png") {
          convert_svg_to_png(image_string);
        } else {
          var pom = document.createElement('a');
          pom.setAttribute('download', 'image.svg');
          pom.setAttribute('href', image_string);
          $("body").append(pom);
          pom.click();
          pom.remove();
        }
      
      }
      
      
      
    

});