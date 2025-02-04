import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonParseException
import com.fasterxml.jackson.databind.JsonMappingException
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper

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
        } catch (IOException e) {
            e.printStackTrace()
        }
    }

    static streamArray(String filePath, ObjectMapper mapper, Closure closure) {
        /*
        Stream the data, one jsonNode at a time from a Json file containing a JsonArray.
        It requires less memory to store maps as JsonNodes than as Maps,
        especially when the maps contain nested maps.
        To use:
        ObjectMapper mapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
        JsonReader.streamArray(filePath, mapper) { JsonNode node ->
            md5 = node.get("md5").asText()
        }
         */
        try {
            JsonFactory factory = mapper.getFactory()
            JsonParser parser = factory.createParser(new File(filePath))
            if (parser.nextToken() == JsonToken.START_ARRAY) {
                while (parser.nextToken() != JsonToken.END_ARRAY) {
                    JsonNode node = mapper.readTree(parser)
                    closure.call(node)
                }
            }
            parser.close()
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $filePath", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath", e)
        }
    }

    static Map<String, JsonNode> load(String filePath, ObjectMapper mapper) {
        /*
        Load all contents of a json file into the memory, storing internal nodes as JsonNodes
        to reduce memory requirements compared to storing the data in memory as maps.
        If can handel json files keyed by a string and valued with a Map.
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
            throw new Exception("File not found: $filePath", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath", e)
        }
    }

    static Map<String, List<JsonNode>> loadMapsWithArrays(String filePath, ObjectMapper mapper) {
        /*
        Load all contents of a JSON file (which is keyed by strings and valued by arrays of JSON objects) into memory,
        storing internal nodes as JsonNodes to reduce memory requirements compared to storing the data as maps.
        Each key in the JSON object maps to an array of JSON objects.
        To use:
        ObjectMapper mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
        def matchesMap = JsonReader.loadMapWithArrayValues(filePath, mapper)
        def arrayNodes = matchesMap.get("EE3D0CF36BFC1CD0F54575DBD337E3AD")
        */
        try {
            File file = new File(filePath)
            JsonNode rootNode = mapper.readTree(file)
            if (!rootNode.isObject()) {
                throw new JsonMappingException("Expected JSON object at root level")
            }

            Map<String, List<JsonNode>> dataMap = new HashMap<>()
            rootNode.fields().each { entry ->
                // If the value is an array, we store the list of JSON nodes
                if (entry.value.isArray()) {
                    List<JsonNode> nodeList = new ArrayList<>()
                    entry.value.each { JsonNode node ->
                        nodeList.add(node)
                    }
                    dataMap.put(entry.key, nodeList)
                } else {
                    // In case the value is not an array, we store it as a single element
                    dataMap.put(entry.key, [entry.value])
                }
            }

            return dataMap
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $filePath", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $filePath", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $filePath", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $filePath", e)
        }
    }
}
