import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.JsonEncoding
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.databind.SerializationFeature

class JsonProcessor {
    private ObjectMapper mapper

    JsonProcessor() {
        this.mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    }

    JsonParser createParser(String jsonFilePath) {
        try {
            JsonParser parser = mapper.getFactory().createParser(new File(jsonFilePath))
            assert parser.nextToken() in [JsonToken.START_OBJECT, JsonToken.START_ARRAY]  // ensure it starts from the beginning of the JSON object/array
            return parser
        } catch (Exception e) {
            throw new Exception("Error creating JsonParser for file $jsonFilePath: ${e}")
        }
    }

    Map jsonToMap(JsonParser parser) {
        try {
            return mapper.readValue(parser, Map)
        } catch (Exception e) {
            throw new Exception("Error reading JSON as Map: ${e}")
        }
    }

    Map jsonToMap(File file) {
        try {
            return mapper.readValue(file, Map)
        } catch (Exception e) {
            throw new Exception("Error reading JSON file as Map: ${e}")
        }
    }

    List jsonToList(JsonParser parser) {
        try {
            return mapper.readValue(parser, List)
        } catch (Exception e) {
            throw new Exception("Error reading JSON as List: ${e}")
        }
    }

    JsonGenerator createGenerator(String outputFilePath) {
        try {
            JsonGenerator generator = mapper.getFactory().createGenerator(new File(outputFilePath), JsonEncoding.UTF8)
            return generator
        } catch (Exception e) {
            throw new Exception("Error creating JsonGenerator for file $outputFilePath: ${e}")
        }
    }

    void write(JsonGenerator generator, Object values) {
        try {
            mapper.writeValue(generator, values)
        } catch (Exception e) {
            throw new Exception("Error writing JSON: ${e}")
        }
    }

    void write(String outputFilePath, values) {
        try {
            mapper.writeValue(new File(outputFilePath), values)
        } catch (Exception e) {
            throw new Exception("Error writing JSON for file $outputFilePath: ${e}")
        }
    }
}