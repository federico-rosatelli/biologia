<script>
export default {
  name: 'NotFound',
};
</script>
<template>
    <div class="not-found">
      <h1 class="title">404 - Page Not Found</h1>
      <p class="message">The page you are looking for does not exist.</p>
    </div>
</template>
  
<style scoped>
    .not-found {
        text-align: center;
    }
    .title {
    font-size: 36px; /* Adjust the font size as desired */
    }

    .message {
    font-size: 18px; /* Adjust the font size as desired */
    }
</style>