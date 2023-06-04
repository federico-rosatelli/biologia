<script>
export default {
	props: ['loading']
}
</script>

<template>
	<div v-if="loading" class="loader">
	</div>
</template>

<style>
    .loader {
        width: 70px;
        height: 70px;
        left: 50%;
        right: 50%;
        border: 3px dotted #FFF;
        border-style: solid solid dotted dotted;
        border-radius: 50%;
        box-sizing: border-box;
        animation: rotation 2s linear infinite;
        position: fixed;
    }

    @keyframes rotation {
        0% {
            transform: rotate(0deg);
        }
        100% {
            transform: rotate(360deg);
        }
    }
    .loader::after {
        content: '';  
        box-sizing: border-box;
        position: absolute;
        left: 0;
        right: 0;
        top: 0;
        bottom: 0;
        margin: auto;
        border: 3px dotted #2a4253;
        border-style: solid solid dotted;
        width: 24px;
        height: 24px;
        border-radius: 50%;
        animation: rotationBack 1s linear infinite;
        transform-origin: center center;
    }
        
    @keyframes rotation {
        0% {
            transform: rotate(0deg);
        }
        100% {
            transform: rotate(360deg);
        }
    } 
    @keyframes rotationBack {
        0% {
            transform: rotate(0deg);
        }
        100% {
            transform: rotate(-360deg);
        }
    } 
</style>