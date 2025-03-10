import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonParseException
import com.fasterxml.jackson.databind.JsonMappingException
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.NullNode

class JsonReader {
    static streamJson(String filePath, ObjectMapper mapper, Closure closure) {
        /*
        Stream through the data in a Json file and pass the field names and values to the closure.
        Streams data as JsonNodes to reduce memory requirements.
        To use:
        ObjectMapper mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
        JsonReader.streamJson(filePath, mapper) { key, value ->
            println("Processed: $key -> $value")
        }
         */
        JsonFactory factory = mapper.getFactory()
        try {
            JsonParser parser = factory.createParser(new File(filePath))
            while (!parser.isClosed()) {
                JsonToken token = parser.nextToken()
                if (token == null) {
                    break
                }

                // Process field names and values
                if (token == JsonToken.FIELD_NAME) {
                    String fieldName = parser.getCurrentName()
                    parser.nextToken()  // Move to the value
                    JsonNode value = mapper.readTree(parser)

                    // Call the closure with the key (field name) and value (JsonNode)
                    closure.call(fieldName, value)
                }
            }
            parser.close()
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        } catch (Exception e) {
            throw new Exception("Error parsing JSON file $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        }
    }
}
