# include <packet.h>

size_t singlePacketToAddress(const Packet to_write, bool * position) {

    bool * start;

    start = position;

    // Write the first 8 flag bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.flag & (0b1 << (7 - i));
        position++;

    }

    // Write the 4 control code bits
    for(int i = 0; i < 4; i++) {

        *position = to_write.controlCode & (0b1 << (3 - i));
        position++;

    }

    // Write the 4 data length code bits
    for(int i = 0; i < 4; i++) {

        *position = to_write.dataLength & (0b1 << (3 - i));
        position++;

    }

    // Write the 8 receiver address bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.recieverAddress & (0b1 << (7 - i));
        position++;

    }

    // Write the 8 sender address bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.senderAddress & (0b1 << (7 - i));
        position++;

    }

    // Write the 8 sequence bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.sequenceNumber & (0b1 << (7 - i));
        position++;

    }

    // Write the data bits
    for(int i = 0; i < to_write.dataLength; i++) {

        for(int j = 0; j < 8; j++) {

            *position = to_write.data.at(i) & (0b1 << (7 - i));
            position++;

        }

    }

    // Write the 32 erc bits
    for(int i = 0; i < 32; i++) {

        *position = to_write.erc & (0b1 << (7 - i));
        position++;

    }

    // Write the final 8 flag bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.flag & (0b1 << (7 - i));
        position++;

    }

    // Return total number of bits written
    return std::distance(start, position);

}

size_t singlePacketToVector(const Packet to_write, std::vector<bool> * buffer) {

    std::vector<bool>::iterator start;
    
    start = buffer->end();

    // Write the first 8 flag bits
    for(int i = 0; i < 8; i++) {

        buffer->push_back(to_write.flag & (0b1 << (7 - i)));

    }

    // Write the 4 control code bits
    for(int i = 0; i < 4; i++) {

        buffer->push_back(to_write.controlCode & (0b1 << (3 - i)));

    }

    // Write the 4 data length code bits
    for(int i = 0; i < 4; i++) {

        buffer->push_back(to_write.dataLength & (0b1 << (3 - i)));

    }

    // Write the 8 receiver address bits
    for(int i = 0; i < 8; i++) {

        buffer->push_back(to_write.recieverAddress & (0b1 << (7 - i)));

    }

    // Write the 8 sender address bits
    for(int i = 0; i < 8; i++) {

        buffer->push_back(to_write.senderAddress & (0b1 << (7 - i)));

    }

    // Write the 8 sequence bits
    for(int i = 0; i < 8; i++) {

        buffer->push_back(to_write.sequenceNumber & (0b1 << (7 - i)));

    }

    // Write the data bits
    for(int i = 0; i < to_write.dataLength; i++) {

        for(int j = 0; j < 8; j++) {

            buffer->push_back(to_write.data.at(i) & (0b1 << (7 - i)));

        }

    }

    // Write the 32 erc bits
    for(int i = 0; i < 32; i++) {

        buffer->push_back(to_write.erc & (0b1 << (7 - i)));

    }

    // Write the final 8 flag bits
    for(int i = 0; i < 8; i++) {

        buffer->push_back(to_write.flag & (0b1 << (7 - i)));

    }

    // Return total number of bits written
    return std::distance(start, buffer->end());

}