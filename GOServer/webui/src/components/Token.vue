<script>
export default {
	props: ["token_value"]
}
</script>

<template>
	<div v-if="token_value">
		<div style="text-align: center">
            <h1>Token success {{token_value}}</h1>
		</div>
	</div>
	<div v-if="!token_value"><slot /></div>
</template>

<style></style>
