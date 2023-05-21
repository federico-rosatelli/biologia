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
                    <div v-if="quals.GBQualifierName != 'transcription'">
                        {{ quals.GBQualifierName }}:
                        {{ quals.GBQualifierValue }}
                    </div>
                </ul>
            </ul>
        </div>
    </div>
</template>
