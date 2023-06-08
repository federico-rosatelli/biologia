<script>
export default {
	props: ['data'],
    methods:{
    }
}
</script>

<template>
    <div>
        {{ data.GBSeqSource }}
        {{ data.GBSeqLocus }}
        <div>
            {{ data.GBSeqComment }}
        </div>
        <div v-for="feature in data.GBSeqFeatureTable">
            <ul>
                {{ feature.GBFeatureKey }}
                <ul v-for="quals in feature.GBFeatureQuals">
                    <div v-if="quals.GBQualifierName != 'transcription'" style="display: flex;">
                        {{ quals.GBQualifierName }}:
                        <a v-if="quals.GBQualifierName == 'product'" :href="'/taxon?product='+quals.GBQualifierValue">
                            {{ quals.GBQualifierValue }}
                        </a>
                        <div v-else>
                            {{ quals.GBQualifierValue }}
                        </div>
                    </div>
                </ul>
            </ul>
        </div>
    </div>
</template>
