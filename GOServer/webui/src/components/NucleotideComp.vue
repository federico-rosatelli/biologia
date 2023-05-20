<script>
export default {
	props: ['data'],
    methods:{
    }
}
</script>

<template>
	<div>
        {{ JSON.stringify(data) }}
    </div>
</template>
