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

    static Map<String, JsonNode> load(String filePath, ObjectMapper mapper) {
        /*
        Load all contents of a json file into the memory, storing internal nodes as JsonNodes
        to reduce memory requirements compared to storing the data in memory as maps.
        It can handel json files keyed by a string and valued with a Map.
        To use:
        ObjectMapper mapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
        def matchesMap = JsonReader.load(filePath, mapper)
        seqMatches = matchesMap.get("seqId_1")
        */
        try {
            File file = new File(filePath)
            JsonNode rootNode = mapper.readTree(file)
            if (!rootNode.isObject()) {
                throw new JsonMappingException("Expected JSON object at root level")
            }

            Map<String, JsonNode> dataMap = new HashMap<>()
            rootNode.fields().each { entry ->
                dataMap.put(entry.key, entry.value)
            }

            return dataMap
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

    static String asString(node) {
        // Convert a textNode to a Groovy string or null without including terminal quotation marks
        if (node instanceof NullNode) {
            return null
        }
        def str = node.asText().replace('\"','')
        return str == "null" ? null : str
    }
}
