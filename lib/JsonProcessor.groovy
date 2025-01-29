import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.JsonEncoding
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.databind.SerializationFeature

class JsonProcessor {
    ObjectMapper mapper
    /*
    Use a single mapper object for all readers and writers to reduce memory requirements.
    Use JsonParser instead of ObjectReader which may perform more internal buffering, thus
    increasing memory requirements.
    */
    JsonProcessor() {
        // Build the main factory class for Jackson
        this.mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    }

    static createParser(String jsonFilePath) {
        // Construct a reader object from the main factory class
        try {
            JsonParser parser = mapper.getFactory().createParser(new File(jsonFilePath))
            assert parser.nextToken() in [JsonToken.START_OBJECT, JsonToken.START_ARRAY]  // ensure it starts from the beginning of the JSON object/array
            return parser
        } catch (Exception e) {
            throw new Exception("Error creating JsonParser for file $jsonFilePath: ${e}")
        }
    }

    static createGenerator(String outputFilePath) {
        // Construct a writer object from the main factory class
        try {
            JsonGenerator generator = mapper.getFactory().createGenerator(new File(outputFilePath), JsonEncoding.UTF8)
            return generator
        } catch (Exception e) {
            throw new Exception("Error creating JsonGenerator for file $outputFilePath: ${e}")
        }
    }

    static jsonToMap(JsonParser parser) {
        try {
            def mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
            return mapper.readValue(parser, Map)
        } catch (Exception e) {
            throw new Exception("Error reading JSON as Map: ${e}")
        }
    }

    static jsonFileToMap(File file) {
        try {
            def mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
            return mapper.readValue(file, Map)
        } catch (Exception e) {
            throw new Exception("Error reading JSON file as Map: ${e}")
        }
    }

    static jsonToList(JsonParser parser) {
        try {
            def mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
            return mapper.readValue(parser, List)
        } catch (Exception e) {
            throw new Exception("Error reading JSON as List: ${e}")
        }
    }

    static generatorWrite(JsonGenerator generator, Object values) {
        try {
            def mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
            mapper.writeValue(generator, values)
        } catch (Exception e) {
            throw new Exception("Error writing JSON: ${e}")
        }
    }

    static write(String outputFilePath, values) {
        try {
            def mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
            mapper.writeValue(new File(outputFilePath), values)
        } catch (Exception e) {
            throw new Exception("Error writing JSON for file $outputFilePath: ${e}")
        }
    }
}