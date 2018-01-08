// code included inside $(document).ready() will only run once the page is ready for JavaScript code to execute
$(document).ready(function() {

  // initialize a counter
  var n = 0;
    // create a click handler which listens for a click on the element with id equal to RStudio
  $("#sis").on("click", function(){

    // increment the counter each time we click on the Rstudio logo
    console.log("meming pretty hard rn");

    // send message to Shiny
    Shiny.onInputChange("meme_dankness", 
        d.map(function(e) {return e.key_var;}));
    });
});
