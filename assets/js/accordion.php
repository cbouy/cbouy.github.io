<script type="text/javascript">
$('.collapse').on('shown.bs.collapse', function(e) {
  var $card = $(this).closest('.card');
  $('html,body').animate({
    scrollTop: $card.offset().top - 55
  }, 100);
});
$(document).ready(function () {
  var hash = location.hash;
  if (hash) {
    var buttonSelected = $(hash).find('button.collapsed');
    if ($(buttonSelected).length) {
      $(buttonSelected).click();
    } else {
      $('html,body').animate({
        scrollTop: $(hash).closest('.card').offset().top - 55
      }, 0);
    }
  }
});
</script>
